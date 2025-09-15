"""
Anomaly Detector for LNB-MDT
============================

This module provides anomaly detection capabilities for molecular dynamics trajectories.
"""

import numpy as np
import pandas as pd
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor
from sklearn.covariance import EllipticEnvelope
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Any, Optional
import logging

class AnomalyDetector:
    """
    Anomaly detection for molecular dynamics data.
    
    This class implements various anomaly detection algorithms to identify
    unusual patterns in molecular dynamics trajectories.
    """
    
    def __init__(self, method: str = 'isolation_forest', **kwargs):
        """
        Initialize the anomaly detector.
        
        Parameters:
        -----------
        method : str
            Anomaly detection method ('isolation_forest', 'lof', 'elliptic_envelope')
        **kwargs : dict
            Additional parameters for the chosen method
        """
        self.method = method
        self.kwargs = kwargs
        self.model = None
        self.scaler = StandardScaler()
        self.is_fitted = False
        
        # Setup logging
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)
        
        # Initialize model based on method
        self._initialize_model()
    
    def _initialize_model(self):
        """Initialize the anomaly detection model."""
        if self.method == 'isolation_forest':
            self.model = IsolationForest(
                contamination=self.kwargs.get('contamination', 0.1),
                random_state=self.kwargs.get('random_state', 42),
                **{k: v for k, v in self.kwargs.items() 
                   if k not in ['contamination', 'random_state']}
            )
        elif self.method == 'lof':
            self.model = LocalOutlierFactor(
                contamination=self.kwargs.get('contamination', 0.1),
                n_neighbors=self.kwargs.get('n_neighbors', 20),
                **{k: v for k, v in self.kwargs.items() 
                   if k not in ['contamination', 'n_neighbors']}
            )
        elif self.method == 'elliptic_envelope':
            self.model = EllipticEnvelope(
                contamination=self.kwargs.get('contamination', 0.1),
                random_state=self.kwargs.get('random_state', 42),
                **{k: v for k, v in self.kwargs.items() 
                   if k not in ['contamination', 'random_state']}
            )
        else:
            raise ValueError(f"Unsupported method: {self.method}")
    
    def fit(self, data: np.ndarray) -> 'AnomalyDetector':
        """
        Fit the anomaly detection model.
        
        Parameters:
        -----------
        data : np.ndarray
            Training data of shape (n_samples, n_features)
        
        Returns:
        --------
        self : AnomalyDetector
        """
        self.logger.info(f"Fitting {self.method} model...")
        
        # Scale the data
        data_scaled = self.scaler.fit_transform(data)
        
        # Fit the model
        if self.method == 'lof':
            # LOF doesn't have a separate fit method
            self.model.fit_predict(data_scaled)
        else:
            self.model.fit(data_scaled)
        
        self.is_fitted = True
        self.logger.info("Model fitting completed")
        
        return self
    
    def predict(self, data: np.ndarray) -> np.ndarray:
        """
        Predict anomalies in the data.
        
        Parameters:
        -----------
        data : np.ndarray
            Data to predict anomalies for
        
        Returns:
        --------
        np.ndarray : Anomaly predictions (-1 for anomalies, 1 for normal)
        """
        if not self.is_fitted:
            raise ValueError("Model must be fitted before prediction")
        
        # Scale the data
        data_scaled = self.scaler.transform(data)
        
        # Make predictions
        if self.method == 'lof':
            predictions = self.model.fit_predict(data_scaled)
        else:
            predictions = self.model.predict(data_scaled)
        
        return predictions
    
    def predict_proba(self, data: np.ndarray) -> np.ndarray:
        """
        Predict anomaly probabilities.
        
        Parameters:
        -----------
        data : np.ndarray
            Data to predict probabilities for
        
        Returns:
        --------
        np.ndarray : Anomaly probabilities
        """
        if not self.is_fitted:
            raise ValueError("Model must be fitted before prediction")
        
        # Scale the data
        data_scaled = self.scaler.transform(data)
        
        # Get decision function scores
        if hasattr(self.model, 'decision_function'):
            scores = self.model.decision_function(data_scaled)
        elif hasattr(self.model, 'score_samples'):
            scores = self.model.score_samples(data_scaled)
        else:
            raise ValueError(f"Model {self.method} doesn't support probability prediction")
        
        # Convert scores to probabilities (higher score = lower anomaly probability)
        probabilities = 1 / (1 + np.exp(-scores))
        
        return probabilities
    
    def detect_anomalies(self, 
                        data: np.ndarray, 
                        threshold: float = 0.5) -> Dict[str, Any]:
        """
        Detect anomalies in the data.
        
        Parameters:
        -----------
        data : np.ndarray
            Data to analyze
        threshold : float
            Probability threshold for anomaly detection
        
        Returns:
        --------
        dict : Detection results
        """
        predictions = self.predict(data)
        probabilities = self.predict_proba(data)
        
        # Identify anomalies
        anomaly_indices = np.where(predictions == -1)[0]
        normal_indices = np.where(predictions == 1)[0]
        
        # Calculate statistics
        n_total = len(data)
        n_anomalies = len(anomaly_indices)
        anomaly_ratio = n_anomalies / n_total
        
        results = {
            'predictions': predictions,
            'probabilities': probabilities,
            'anomaly_indices': anomaly_indices,
            'normal_indices': normal_indices,
            'n_anomalies': n_anomalies,
            'n_normal': len(normal_indices),
            'anomaly_ratio': anomaly_ratio,
            'threshold': threshold
        }
        
        self.logger.info(f"Detected {n_anomalies} anomalies out of {n_total} samples ({anomaly_ratio:.2%})")
        
        return results


class MDAnomalyDetector(AnomalyDetector):
    """
    Specialized anomaly detector for molecular dynamics data.
    """
    
    def __init__(self, **kwargs):
        """Initialize MD-specific anomaly detector."""
        super().__init__(**kwargs)
        self.feature_names = []
    
    def extract_features(self, 
                        trajectory_data: Dict[str, np.ndarray],
                        feature_types: List[str] = None) -> np.ndarray:
        """
        Extract features from molecular dynamics trajectory data.
        
        Parameters:
        -----------
        trajectory_data : dict
            Dictionary containing trajectory data (e.g., positions, velocities, etc.)
        feature_types : list
            Types of features to extract
        
        Returns:
        --------
        np.ndarray : Extracted features
        """
        if feature_types is None:
            feature_types = ['position_stats', 'velocity_stats', 'distance_stats']
        
        features = []
        self.feature_names = []
        
        for feature_type in feature_types:
            if feature_type == 'position_stats':
                pos_features = self._extract_position_features(trajectory_data)
                features.append(pos_features)
                self.feature_names.extend([f'pos_{i}' for i in range(pos_features.shape[1])])
            
            elif feature_type == 'velocity_stats':
                vel_features = self._extract_velocity_features(trajectory_data)
                features.append(vel_features)
                self.feature_names.extend([f'vel_{i}' for i in range(vel_features.shape[1])])
            
            elif feature_type == 'distance_stats':
                dist_features = self._extract_distance_features(trajectory_data)
                features.append(dist_features)
                self.feature_names.extend([f'dist_{i}' for i in range(dist_features.shape[1])])
        
        return np.hstack(features)
    
    def _extract_position_features(self, trajectory_data: Dict[str, np.ndarray]) -> np.ndarray:
        """Extract position-based features."""
        if 'positions' not in trajectory_data:
            return np.array([])
        
        positions = trajectory_data['positions']  # Shape: (n_frames, n_atoms, 3)
        
        # Calculate statistics for each frame
        features = []
        for frame in positions:
            # Mean position
            mean_pos = np.mean(frame, axis=0)
            
            # Standard deviation of positions
            std_pos = np.std(frame, axis=0)
            
            # Range of positions
            pos_range = np.ptp(frame, axis=0)
            
            # Combine features
            frame_features = np.concatenate([mean_pos, std_pos, pos_range])
            features.append(frame_features)
        
        return np.array(features)
    
    def _extract_velocity_features(self, trajectory_data: Dict[str, np.ndarray]) -> np.ndarray:
        """Extract velocity-based features."""
        if 'velocities' not in trajectory_data:
            return np.array([])
        
        velocities = trajectory_data['velocities']  # Shape: (n_frames, n_atoms, 3)
        
        # Calculate statistics for each frame
        features = []
        for frame in velocities:
            # Mean velocity
            mean_vel = np.mean(frame, axis=0)
            
            # Standard deviation of velocities
            std_vel = np.std(frame, axis=0)
            
            # Kinetic energy
            kinetic_energy = 0.5 * np.sum(frame**2, axis=1)
            mean_ke = np.mean(kinetic_energy)
            std_ke = np.std(kinetic_energy)
            
            # Combine features
            frame_features = np.concatenate([mean_vel, std_vel, [mean_ke, std_ke]])
            features.append(frame_features)
        
        return np.array(features)
    
    def _extract_distance_features(self, trajectory_data: Dict[str, np.ndarray]) -> np.ndarray:
        """Extract distance-based features."""
        if 'positions' not in trajectory_data:
            return np.array([])
        
        positions = trajectory_data['positions']  # Shape: (n_frames, n_atoms, 3)
        
        # Calculate statistics for each frame
        features = []
        for frame in positions:
            # Calculate pairwise distances
            distances = []
            for i in range(len(frame)):
                for j in range(i+1, len(frame)):
                    dist = np.linalg.norm(frame[i] - frame[j])
                    distances.append(dist)
            
            distances = np.array(distances)
            
            # Distance statistics
            mean_dist = np.mean(distances)
            std_dist = np.std(distances)
            min_dist = np.min(distances)
            max_dist = np.max(distances)
            
            # Combine features
            frame_features = np.array([mean_dist, std_dist, min_dist, max_dist])
            features.append(frame_features)
        
        return np.array(features)
    
    def analyze_trajectory(self, 
                          gro_file: str,
                          xtc_file: str,
                          residues: Dict[str, List[str]],
                          start_frame: int = 0,
                          stop_frame: int = None,
                          step_frame: int = 1) -> Dict[str, Any]:
        """
        Analyze a molecular dynamics trajectory for anomalies.
        
        Parameters:
        -----------
        gro_file : str
            Path to GRO file
        xtc_file : str
            Path to XTC file
        residues : dict
            Residue groups to analyze
        start_frame : int
            Starting frame
        stop_frame : int
            Stopping frame
        step_frame : int
            Frame step size
        
        Returns:
        --------
        dict : Analysis results
        """
        try:
            import MDAnalysis as mda
            
            # Load trajectory
            universe = mda.Universe(gro_file, xtc_file)
            
            # Select atoms based on residues
            selection = universe.select_atoms(' or '.join([
                f'resname {resname} and name {" ".join(names)}'
                for resname, names in residues.items()
            ]))
            
            # Extract trajectory data
            trajectory_data = self._extract_trajectory_data(
                universe, selection, start_frame, stop_frame, step_frame
            )
            
            # Extract features
            features = self.extract_features(trajectory_data)
            
            # Fit model and detect anomalies
            self.fit(features)
            results = self.detect_anomalies(features)
            
            # Add trajectory information
            results['trajectory_info'] = {
                'n_frames': len(features),
                'n_atoms': len(selection),
                'residues': residues,
                'feature_names': self.feature_names
            }
            
            return results
            
        except ImportError:
            self.logger.error("MDAnalysis not available")
            raise
        except Exception as e:
            self.logger.error(f"Error analyzing trajectory: {e}")
            raise
    
    def _extract_trajectory_data(self, 
                                universe: 'mda.Universe',
                                selection: 'mda.AtomGroup',
                                start_frame: int,
                                stop_frame: int,
                                step_frame: int) -> Dict[str, np.ndarray]:
        """Extract data from MDAnalysis trajectory."""
        positions = []
        velocities = []
        
        for ts in universe.trajectory[start_frame:stop_frame:step_frame]:
            positions.append(selection.positions.copy())
            if hasattr(selection, 'velocities') and selection.velocities is not None:
                velocities.append(selection.velocities.copy())
        
        trajectory_data = {'positions': np.array(positions)}
        if velocities:
            trajectory_data['velocities'] = np.array(velocities)
        
        return trajectory_data
    
    def plot_anomalies(self, 
                      results: Dict[str, Any],
                      save_path: str = None) -> None:
        """
        Plot anomaly detection results.
        
        Parameters:
        -----------
        results : dict
            Results from detect_anomalies
        save_path : str
            Path to save the plot
        """
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Plot 1: Anomaly predictions over time
        axes[0, 0].plot(results['predictions'], 'b-', alpha=0.7, label='Predictions')
        axes[0, 0].scatter(results['anomaly_indices'], 
                          results['predictions'][results['anomaly_indices']], 
                          color='red', s=50, label='Anomalies')
        axes[0, 0].set_title('Anomaly Predictions Over Time')
        axes[0, 0].set_xlabel('Frame')
        axes[0, 0].set_ylabel('Prediction (-1: Anomaly, 1: Normal)')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # Plot 2: Anomaly probabilities
        axes[0, 1].plot(results['probabilities'], 'g-', alpha=0.7)
        axes[0, 1].axhline(y=results['threshold'], color='red', linestyle='--', 
                          label=f'Threshold ({results["threshold"]})')
        axes[0, 1].set_title('Anomaly Probabilities')
        axes[0, 1].set_xlabel('Frame')
        axes[0, 1].set_ylabel('Anomaly Probability')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        # Plot 3: Anomaly distribution
        axes[1, 0].hist(results['probabilities'], bins=30, alpha=0.7, color='skyblue')
        axes[1, 0].axvline(x=results['threshold'], color='red', linestyle='--', 
                          label=f'Threshold ({results["threshold"]})')
        axes[1, 0].set_title('Distribution of Anomaly Probabilities')
        axes[1, 0].set_xlabel('Anomaly Probability')
        axes[1, 0].set_ylabel('Frequency')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        # Plot 4: Summary statistics
        summary_text = f"""
        Total Frames: {len(results['predictions'])}
        Anomalies: {results['n_anomalies']}
        Normal: {results['n_normal']}
        Anomaly Ratio: {results['anomaly_ratio']:.2%}
        Method: {self.method}
        """
        axes[1, 1].text(0.1, 0.5, summary_text, transform=axes[1, 1].transAxes,
                       fontsize=12, verticalalignment='center',
                       bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
        axes[1, 1].set_title('Summary Statistics')
        axes[1, 1].axis('off')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            self.logger.info(f"Plot saved to {save_path}")
        
        plt.show()


# Example usage
if __name__ == "__main__":
    # Example: Detect anomalies in trajectory
    detector = MDAnomalyDetector(method='isolation_forest', contamination=0.1)
    
    # Analyze trajectory
    results = detector.analyze_trajectory(
        gro_file="cases/lnb.gro",
        xtc_file="cases/md.xtc",
        residues={'DPPC': ['PO4'], 'CHOL': ['ROH']},
        start_frame=0,
        stop_frame=100,
        step_frame=1
    )
    
    # Plot results
    detector.plot_anomalies(results, save_path="anomaly_detection_results.png")
    
    print(f"Detected {results['n_anomalies']} anomalies")
    print(f"Anomaly ratio: {results['anomaly_ratio']:.2%}")
