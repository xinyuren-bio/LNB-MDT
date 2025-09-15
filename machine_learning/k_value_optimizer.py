"""
K-Value Optimizer for LNB-MDT
=============================

This module provides intelligent k-value optimization for molecular dynamics analysis.
"""

import numpy as np
import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Any, Optional, Callable
import logging
import time
import joblib

class KValueOptimizer:
    """
    Intelligent k-value optimizer for molecular dynamics analysis.
    
    This class uses machine learning to automatically find the optimal k-value
    for different analysis types based on system properties.
    """
    
    def __init__(self, 
                 analysis_type: str,
                 system_properties: Dict[str, float] = None,
                 k_range: Tuple[int, int] = (5, 50),
                 n_trials: int = 20,
                 random_state: int = 42):
        """
        Initialize the k-value optimizer.
        
        Parameters:
        -----------
        analysis_type : str
            Type of analysis ('area', 'curvature', 'height', 'cluster')
        system_properties : dict
            System properties (density, temperature, pressure, etc.)
        k_range : tuple
            Range of k-values to test (min, max)
        n_trials : int
            Number of trials for optimization
        random_state : int
            Random seed for reproducibility
        """
        self.analysis_type = analysis_type
        self.system_properties = system_properties or {}
        self.k_range = k_range
        self.n_trials = n_trials
        self.random_state = random_state
        
        # Initialize Gaussian Process
        kernel = ConstantKernel(1.0) * RBF(length_scale=1.0)
        self.gp = GaussianProcessRegressor(
            kernel=kernel,
            alpha=1e-6,
            normalize_y=True,
            random_state=random_state
        )
        
        # Storage for optimization history
        self.k_values_tested = []
        self.objective_values = []
        self.best_k = None
        self.best_score = float('inf')
        
        # Setup logging
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)
    
    def _get_default_k_bounds(self) -> Dict[str, Tuple[int, int]]:
        """Get default k-value bounds for different analysis types."""
        bounds = {
            'area': (10, 40),        # Voronoi tessellation
            'curvature': (15, 35),   # Surface curvature
            'height': (8, 30),       # Height analysis
            'cluster': (5, 25),      # Clustering
            'pca': (20, 50),         # PCA analysis
            'anisotropy': (12, 35),  # Anisotropy
            'gyration': (10, 30),    # Gyration radius
            'sz': (15, 40)           # Sz order parameter
        }
        return bounds.get(self.analysis_type, self.k_range)
    
    def _extract_system_features(self, gro_file: str, xtc_file: str) -> Dict[str, float]:
        """
        Extract system features that influence optimal k-value.
        
        Parameters:
        -----------
        gro_file : str
            Path to GRO file
        xtc_file : str
            Path to XTC file
        
        Returns:
        --------
        dict : System features
        """
        try:
            import MDAnalysis as mda
            
            # Load system
            universe = mda.Universe(gro_file, xtc_file)
            
            # Extract system properties
            features = {}
            
            # Number of atoms
            features['n_atoms'] = len(universe.atoms)
            
            # System size
            if hasattr(universe, 'dimensions'):
                box = universe.dimensions[:3]
                features['box_volume'] = np.prod(box)
                features['box_length'] = np.mean(box)
            else:
                features['box_volume'] = 1000.0  # Default
                features['box_length'] = 10.0    # Default
            
            # Density
            features['density'] = features['n_atoms'] / features['box_volume']
            
            # Number of frames
            features['n_frames'] = len(universe.trajectory)
            
            # Residue types (if available)
            if hasattr(universe, 'residues'):
                features['n_residues'] = len(universe.residues)
                features['avg_residue_size'] = features['n_atoms'] / features['n_residues']
            else:
                features['n_residues'] = features['n_atoms'] // 10  # Estimate
                features['avg_residue_size'] = 10.0
            
            # Add user-provided properties
            features.update(self.system_properties)
            
            return features
            
        except ImportError:
            self.logger.warning("MDAnalysis not available, using default features")
            return {
                'n_atoms': 1000,
                'density': 0.1,
                'box_length': 10.0,
                'n_residues': 100,
                'avg_residue_size': 10.0
            }
    
    def _calculate_optimal_k_heuristic(self, features: Dict[str, float]) -> int:
        """
        Calculate initial k-value using heuristic rules.
        
        Parameters:
        -----------
        features : dict
            System features
        
        Returns:
        --------
        int : Initial k-value estimate
        """
        # Base k-value based on analysis type
        base_k = {
            'area': 20,
            'curvature': 25,
            'height': 15,
            'cluster': 10,
            'pca': 30,
            'anisotropy': 20,
            'gyration': 15,
            'sz': 25
        }.get(self.analysis_type, 20)
        
        # Adjust based on system density
        density_factor = np.clip(features.get('density', 0.1) / 0.1, 0.5, 2.0)
        
        # Adjust based on system size
        size_factor = np.clip(features.get('n_atoms', 1000) / 1000, 0.5, 2.0)
        
        # Calculate optimal k
        optimal_k = int(base_k * density_factor * size_factor)
        
        # Ensure within bounds
        k_min, k_max = self._get_default_k_bounds()
        optimal_k = np.clip(optimal_k, k_min, k_max)
        
        return optimal_k
    
    def _evaluate_k_value(self, k_value: int, gro_file: str, xtc_file: str, 
                         residues: Dict[str, List[str]]) -> float:
        """
        Evaluate a specific k-value by running the analysis.
        
        Parameters:
        -----------
        k_value : int
            K-value to test
        gro_file : str
            Path to GRO file
        xtc_file : str
            Path to XTC file
        residues : dict
            Residue groups for analysis
        
        Returns:
        --------
        float : Objective value (lower is better)
        """
        try:
            # Import analysis module dynamically
            if self.analysis_type == 'area':
                from analysis.area import Area
                analyzer = Area(gro_file, xtc_file, residues, k_value=k_value)
            elif self.analysis_type == 'curvature':
                from analysis.curvature import Curvature
                analyzer = Curvature(gro_file, xtc_file, residues, k_value=k_value)
            elif self.analysis_type == 'height':
                from analysis.height import Height
                import MDAnalysis as mda
                universe = mda.Universe(gro_file, xtc_file)
                analyzer = Height(universe=universe, residuesGroup=residues, k=k_value)
            elif self.analysis_type == 'cluster':
                from analysis.cluster import Cluster
                analyzer = Cluster(gro_file, xtc_file, residues, cutoff=k_value)
            else:
                raise ValueError(f"Unsupported analysis type: {self.analysis_type}")
            
            # Run analysis and measure performance
            start_time = time.time()
            results = analyzer.run()
            computation_time = time.time() - start_time
            
            # Calculate objective value
            # Consider computation time, result quality, and stability
            objective = self._calculate_objective(k_value, results, computation_time)
            
            return objective
            
        except Exception as e:
            self.logger.warning(f"Error evaluating k={k_value}: {e}")
            return float('inf')  # High penalty for failed evaluations
    
    def _calculate_objective(self, k_value: int, results: Any, computation_time: float) -> float:
        """
        Calculate objective value for optimization.
        
        Parameters:
        -----------
        k_value : int
            K-value used
        results : Any
            Analysis results
        computation_time : float
            Time taken for computation
        
        Returns:
        --------
        float : Objective value (lower is better)
        """
        # Base objective: computation time
        objective = computation_time
        
        # Penalty for very small k (poor accuracy)
        if k_value < 10:
            objective += (10 - k_value) * 0.1
        
        # Penalty for very large k (overfitting)
        if k_value > 40:
            objective += (k_value - 40) * 0.05
        
        # Quality penalty based on result characteristics
        if hasattr(results, '__len__'):
            # Penalty for too few or too many results
            n_results = len(results)
            if n_results < 10:
                objective += (10 - n_results) * 0.01
            elif n_results > 1000:
                objective += (n_results - 1000) * 0.001
        
        return objective
    
    def optimize(self, gro_file: str, xtc_file: str, 
                residues: Dict[str, List[str]]) -> Dict[str, Any]:
        """
        Optimize k-value for the given system and analysis.
        
        Parameters:
        -----------
        gro_file : str
            Path to GRO file
        xtc_file : str
            Path to XTC file
        residues : dict
            Residue groups for analysis
        
        Returns:
        --------
        dict : Optimization results
        """
        self.logger.info(f"Starting k-value optimization for {self.analysis_type} analysis...")
        
        # Extract system features
        features = self._extract_system_features(gro_file, xtc_file)
        self.logger.info(f"System features: {features}")
        
        # Get k-value bounds
        k_min, k_max = self._get_default_k_bounds()
        
        # Calculate initial k-value using heuristic
        initial_k = self._calculate_optimal_k_heuristic(features)
        self.logger.info(f"Initial k-value estimate: {initial_k}")
        
        # Test initial k-value
        initial_score = self._evaluate_k_value(initial_k, gro_file, xtc_file, residues)
        self.k_values_tested.append(initial_k)
        self.objective_values.append(initial_score)
        
        if initial_score < self.best_score:
            self.best_score = initial_score
            self.best_k = initial_k
        
        # Bayesian optimization
        for trial in range(self.n_trials - 1):
            self.logger.info(f"Optimization trial {trial + 1}/{self.n_trials - 1}")
            
            # Fit Gaussian Process
            X = np.array(self.k_values_tested).reshape(-1, 1)
            y = np.array(self.objective_values)
            self.gp.fit(X, y)
            
            # Find next k-value to test
            k_candidates = np.arange(k_min, k_max + 1)
            k_candidates = k_candidates[~np.isin(k_candidates, self.k_values_tested)]
            
            if len(k_candidates) == 0:
                break
            
            # Predict objective values
            k_candidates_2d = k_candidates.reshape(-1, 1)
            mean_pred, std_pred = self.gp.predict(k_candidates_2d, return_std=True)
            
            # Expected Improvement acquisition function
            best_y = min(self.objective_values)
            improvement = best_y - mean_pred
            z = improvement / (std_pred + 1e-8)
            ei = improvement * self._normal_cdf(z) + std_pred * self._normal_pdf(z)
            
            # Select k-value with maximum expected improvement
            next_k = k_candidates[np.argmax(ei)]
            
            # Evaluate next k-value
            score = self._evaluate_k_value(next_k, gro_file, xtc_file, residues)
            
            self.k_values_tested.append(next_k)
            self.objective_values.append(score)
            
            if score < self.best_score:
                self.best_score = score
                self.best_k = next_k
                self.logger.info(f"New best k-value: {next_k} (score: {score:.4f})")
        
        return {
            'best_k': self.best_k,
            'best_score': self.best_score,
            'k_values_tested': self.k_values_tested,
            'objective_values': self.objective_values,
            'system_features': features,
            'analysis_type': self.analysis_type
        }
    
    def _normal_cdf(self, x):
        """Normal cumulative distribution function."""
        return 0.5 * (1 + np.sign(x) * np.sqrt(1 - np.exp(-2 * x**2 / np.pi)))
    
    def _normal_pdf(self, x):
        """Normal probability density function."""
        return np.exp(-0.5 * x**2) / np.sqrt(2 * np.pi)
    
    def plot_optimization_results(self, results: Dict[str, Any], save_path: str = None):
        """
        Plot optimization results.
        
        Parameters:
        -----------
        results : dict
            Optimization results
        save_path : str
            Path to save the plot
        """
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Plot 1: Objective values vs k-values
        axes[0, 0].scatter(results['k_values_tested'], results['objective_values'], 
                          alpha=0.7, color='blue', s=50)
        axes[0, 0].axvline(x=results['best_k'], color='red', linestyle='--', 
                          label=f'Best k = {results["best_k"]}')
        axes[0, 0].set_xlabel('K-Value')
        axes[0, 0].set_ylabel('Objective Value')
        axes[0, 0].set_title('K-Value Optimization Results')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # Plot 2: Optimization progress
        axes[0, 1].plot(range(len(results['objective_values'])), 
                       np.minimum.accumulate(results['objective_values']), 
                       'g-', linewidth=2)
        axes[0, 1].set_xlabel('Trial Number')
        axes[0, 1].set_ylabel('Best Objective Value')
        axes[0, 1].set_title('Optimization Progress')
        axes[0, 1].grid(True, alpha=0.3)
        
        # Plot 3: K-value distribution
        axes[1, 0].hist(results['k_values_tested'], bins=10, alpha=0.7, color='skyblue')
        axes[1, 0].axvline(x=results['best_k'], color='red', linestyle='--', 
                          label=f'Best k = {results["best_k"]}')
        axes[1, 0].set_xlabel('K-Value')
        axes[1, 0].set_ylabel('Frequency')
        axes[1, 0].set_title('K-Value Distribution')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        # Plot 4: System features summary
        features = results['system_features']
        feature_names = list(features.keys())[:5]  # Show top 5 features
        feature_values = [features[name] for name in feature_names]
        
        axes[1, 1].barh(range(len(feature_names)), feature_values)
        axes[1, 1].set_yticks(range(len(feature_names)))
        axes[1, 1].set_yticklabels(feature_names)
        axes[1, 1].set_xlabel('Value')
        axes[1, 1].set_title('System Features')
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            self.logger.info(f"Plot saved to {save_path}")
        
        plt.show()
    
    def get_recommendations(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """
        Get k-value recommendations based on optimization results.
        
        Parameters:
        -----------
        results : dict
            Optimization results
        
        Returns:
        --------
        dict : Recommendations
        """
        features = results['system_features']
        best_k = results['best_k']
        
        recommendations = {
            'optimal_k': best_k,
            'confidence': 'high' if len(results['k_values_tested']) > 10 else 'medium',
            'reasoning': [],
            'alternative_k_values': []
        }
        
        # Generate reasoning
        if features.get('density', 0) > 0.15:
            recommendations['reasoning'].append("High system density requires larger k-value")
        elif features.get('density', 0) < 0.05:
            recommendations['reasoning'].append("Low system density allows smaller k-value")
        
        if features.get('n_atoms', 0) > 2000:
            recommendations['reasoning'].append("Large system benefits from larger k-value")
        elif features.get('n_atoms', 0) < 500:
            recommendations['reasoning'].append("Small system works well with smaller k-value")
        
        # Suggest alternative k-values
        k_min, k_max = self._get_default_k_bounds()
        recommendations['alternative_k_values'] = [
            max(k_min, best_k - 2),
            min(k_max, best_k + 2),
            max(k_min, best_k - 5),
            min(k_max, best_k + 5)
        ]
        
        return recommendations


class SmartKValueSelector:
    """
    Smart k-value selector that learns from previous optimizations.
    """
    
    def __init__(self):
        """Initialize the smart selector."""
        self.optimization_history = []
        self.model = None
        self.scaler = StandardScaler()
        self.is_trained = False
    
    def add_optimization_result(self, result: Dict[str, Any]):
        """Add optimization result to history."""
        self.optimization_history.append(result)
    
    def train_model(self):
        """Train a model to predict optimal k-values."""
        if len(self.optimization_history) < 5:
            self.logger.warning("Need at least 5 optimization results to train model")
            return
        
        # Prepare training data
        X = []
        y = []
        
        for result in self.optimization_history:
            features = result['system_features']
            optimal_k = result['best_k']
            
            # Create feature vector
            feature_vector = [
                features.get('n_atoms', 1000),
                features.get('density', 0.1),
                features.get('box_length', 10.0),
                features.get('n_residues', 100),
                features.get('avg_residue_size', 10.0)
            ]
            
            X.append(feature_vector)
            y.append(optimal_k)
        
        X = np.array(X)
        y = np.array(y)
        
        # Scale features
        X_scaled = self.scaler.fit_transform(X)
        
        # Train model (simple linear regression for now)
        from sklearn.linear_model import LinearRegression
        self.model = LinearRegression()
        self.model.fit(X_scaled, y)
        
        self.is_trained = True
        self.logger.info("Smart k-value selector trained successfully")
    
    def predict_optimal_k(self, features: Dict[str, float]) -> int:
        """
        Predict optimal k-value for new system.
        
        Parameters:
        -----------
        features : dict
            System features
        
        Returns:
        --------
        int : Predicted optimal k-value
        """
        if not self.is_trained:
            raise ValueError("Model must be trained before prediction")
        
        # Create feature vector
        feature_vector = [
            features.get('n_atoms', 1000),
            features.get('density', 0.1),
            features.get('box_length', 10.0),
            features.get('n_residues', 100),
            features.get('avg_residue_size', 10.0)
        ]
        
        # Scale features
        X_scaled = self.scaler.transform([feature_vector])
        
        # Predict
        predicted_k = int(self.model.predict(X_scaled)[0])
        
        return max(5, min(50, predicted_k))  # Ensure within reasonable bounds


# Example usage
if __name__ == "__main__":
    # Example: Optimize k-value for area analysis
    optimizer = KValueOptimizer('area', n_trials=15)
    
    # Run optimization
    results = optimizer.optimize(
        gro_file="cases/lnb.gro",
        xtc_file="cases/md.xtc",
        residues={'DPPC': ['PO4'], 'CHOL': ['ROH']}
    )
    
    print(f"Optimal k-value: {results['best_k']}")
    print(f"Best score: {results['best_score']:.4f}")
    
    # Plot results
    optimizer.plot_optimization_results(results, save_path="k_value_optimization.png")
    
    # Get recommendations
    recommendations = optimizer.get_recommendations(results)
    print(f"Recommendations: {recommendations}")
