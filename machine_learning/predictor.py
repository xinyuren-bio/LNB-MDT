"""
Property Predictor for LNB-MDT
==============================

This module provides predictive modeling capabilities for molecular properties.
"""

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.linear_model import LinearRegression, Ridge, Lasso
from sklearn.svm import SVR
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
import joblib
from typing import Dict, List, Tuple, Any, Optional, Union
import logging
import matplotlib.pyplot as plt
import seaborn as sns

class PropertyPredictor:
    """
    Property predictor for molecular dynamics data.
    
    This class implements various machine learning models to predict
    molecular properties from trajectory data.
    """
    
    def __init__(self, 
                 model_type: str = 'random_forest',
                 target_property: str = 'diffusion_coefficient',
                 **kwargs):
        """
        Initialize the property predictor.
        
        Parameters:
        -----------
        model_type : str
            Type of model ('random_forest', 'gradient_boosting', 'linear', 'svr', 'neural_network')
        target_property : str
            Property to predict
        **kwargs : dict
            Additional parameters for the model
        """
        self.model_type = model_type
        self.target_property = target_property
        self.kwargs = kwargs
        self.model = None
        self.scaler = StandardScaler()
        self.is_fitted = False
        self.feature_names = []
        
        # Setup logging
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)
        
        # Initialize model
        self._initialize_model()
    
    def _initialize_model(self):
        """Initialize the machine learning model."""
        if self.model_type == 'random_forest':
            self.model = RandomForestRegressor(
                n_estimators=self.kwargs.get('n_estimators', 100),
                max_depth=self.kwargs.get('max_depth', None),
                random_state=self.kwargs.get('random_state', 42),
                **{k: v for k, v in self.kwargs.items() 
                   if k not in ['n_estimators', 'max_depth', 'random_state']}
            )
        elif self.model_type == 'gradient_boosting':
            self.model = GradientBoostingRegressor(
                n_estimators=self.kwargs.get('n_estimators', 100),
                learning_rate=self.kwargs.get('learning_rate', 0.1),
                max_depth=self.kwargs.get('max_depth', 3),
                random_state=self.kwargs.get('random_state', 42),
                **{k: v for k, v in self.kwargs.items() 
                   if k not in ['n_estimators', 'learning_rate', 'max_depth', 'random_state']}
            )
        elif self.model_type == 'linear':
            self.model = LinearRegression(**self.kwargs)
        elif self.model_type == 'ridge':
            self.model = Ridge(
                alpha=self.kwargs.get('alpha', 1.0),
                random_state=self.kwargs.get('random_state', 42),
                **{k: v for k, v in self.kwargs.items() 
                   if k not in ['alpha', 'random_state']}
            )
        elif self.model_type == 'lasso':
            self.model = Lasso(
                alpha=self.kwargs.get('alpha', 1.0),
                random_state=self.kwargs.get('random_state', 42),
                **{k: v for k, v in self.kwargs.items() 
                   if k not in ['alpha', 'random_state']}
            )
        elif self.model_type == 'svr':
            self.model = SVR(
                kernel=self.kwargs.get('kernel', 'rbf'),
                C=self.kwargs.get('C', 1.0),
                gamma=self.kwargs.get('gamma', 'scale'),
                **{k: v for k, v in self.kwargs.items() 
                   if k not in ['kernel', 'C', 'gamma']}
            )
        elif self.model_type == 'neural_network':
            self.model = MLPRegressor(
                hidden_layer_sizes=self.kwargs.get('hidden_layer_sizes', (100, 50)),
                activation=self.kwargs.get('activation', 'relu'),
                solver=self.kwargs.get('solver', 'adam'),
                alpha=self.kwargs.get('alpha', 0.0001),
                random_state=self.kwargs.get('random_state', 42),
                **{k: v for k, v in self.kwargs.items() 
                   if k not in ['hidden_layer_sizes', 'activation', 'solver', 'alpha', 'random_state']}
            )
        else:
            raise ValueError(f"Unsupported model type: {self.model_type}")
    
    def extract_features(self, 
                        trajectory_data: Dict[str, np.ndarray],
                        feature_types: List[str] = None) -> np.ndarray:
        """
        Extract features from trajectory data.
        
        Parameters:
        -----------
        trajectory_data : dict
            Dictionary containing trajectory data
        feature_types : list
            Types of features to extract
        
        Returns:
        --------
        np.ndarray : Extracted features
        """
        if feature_types is None:
            feature_types = ['structural', 'dynamical', 'thermodynamic']
        
        features = []
        self.feature_names = []
        
        for feature_type in feature_types:
            if feature_type == 'structural':
                struct_features = self._extract_structural_features(trajectory_data)
                features.append(struct_features)
                self.feature_names.extend([f'struct_{i}' for i in range(struct_features.shape[1])])
            
            elif feature_type == 'dynamical':
                dyn_features = self._extract_dynamical_features(trajectory_data)
                features.append(dyn_features)
                self.feature_names.extend([f'dyn_{i}' for i in range(dyn_features.shape[1])])
            
            elif feature_type == 'thermodynamic':
                thermo_features = self._extract_thermodynamic_features(trajectory_data)
                features.append(thermo_features)
                self.feature_names.extend([f'thermo_{i}' for i in range(thermo_features.shape[1])])
        
        return np.hstack(features)
    
    def _extract_structural_features(self, trajectory_data: Dict[str, np.ndarray]) -> np.ndarray:
        """Extract structural features."""
        if 'positions' not in trajectory_data:
            return np.array([])
        
        positions = trajectory_data['positions']
        
        # Calculate structural features for each frame
        features = []
        for frame in positions:
            # Radius of gyration
            center_of_mass = np.mean(frame, axis=0)
            rg_squared = np.mean(np.sum((frame - center_of_mass)**2, axis=1))
            rg = np.sqrt(rg_squared)
            
            # End-to-end distance
            if len(frame) >= 2:
                end_to_end = np.linalg.norm(frame[-1] - frame[0])
            else:
                end_to_end = 0.0
            
            # Asphericity
            if len(frame) >= 3:
                # Calculate gyration tensor
                gyration_tensor = np.zeros((3, 3))
                for atom_pos in frame:
                    r = atom_pos - center_of_mass
                    gyration_tensor += np.outer(r, r)
                gyration_tensor /= len(frame)
                
                # Eigenvalues
                eigenvals = np.linalg.eigvals(gyration_tensor)
                eigenvals = np.sort(eigenvals.real)[::-1]  # Sort in descending order
                
                # Asphericity
                asphericity = eigenvals[0] - 0.5 * (eigenvals[1] + eigenvals[2])
            else:
                asphericity = 0.0
            
            frame_features = np.array([rg, end_to_end, asphericity])
            features.append(frame_features)
        
        return np.array(features)
    
    def _extract_dynamical_features(self, trajectory_data: Dict[str, np.ndarray]) -> np.ndarray:
        """Extract dynamical features."""
        if 'velocities' not in trajectory_data:
            return np.array([])
        
        velocities = trajectory_data['velocities']
        
        # Calculate dynamical features for each frame
        features = []
        for frame in velocities:
            # Mean velocity
            mean_vel = np.mean(frame, axis=0)
            mean_vel_magnitude = np.linalg.norm(mean_vel)
            
            # Velocity fluctuations
            vel_fluctuations = frame - mean_vel
            vel_fluctuation_magnitude = np.mean(np.linalg.norm(vel_fluctuations, axis=1))
            
            # Kinetic energy
            kinetic_energy = 0.5 * np.sum(frame**2, axis=1)
            mean_ke = np.mean(kinetic_energy)
            ke_std = np.std(kinetic_energy)
            
            frame_features = np.array([mean_vel_magnitude, vel_fluctuation_magnitude, mean_ke, ke_std])
            features.append(frame_features)
        
        return np.array(features)
    
    def _extract_thermodynamic_features(self, trajectory_data: Dict[str, np.ndarray]) -> np.ndarray:
        """Extract thermodynamic features."""
        if 'positions' not in trajectory_data:
            return np.array([])
        
        positions = trajectory_data['positions']
        
        # Calculate thermodynamic features for each frame
        features = []
        for frame in positions:
            # Volume (approximate using convex hull)
            try:
                from scipy.spatial import ConvexHull
                hull = ConvexHull(frame)
                volume = hull.volume
            except:
                volume = 0.0
            
            # Surface area
            try:
                surface_area = hull.area
            except:
                surface_area = 0.0
            
            # Density
            density = len(frame) / (volume + 1e-8)
            
            frame_features = np.array([volume, surface_area, density])
            features.append(frame_features)
        
        return np.array(features)
    
    def fit(self, X: np.ndarray, y: np.ndarray, 
            test_size: float = 0.2, 
            random_state: int = 42) -> Dict[str, Any]:
        """
        Fit the model to the data.
        
        Parameters:
        -----------
        X : np.ndarray
            Feature matrix
        y : np.ndarray
            Target values
        test_size : float
            Fraction of data to use for testing
        random_state : int
            Random seed for reproducibility
        
        Returns:
        --------
        dict : Training results
        """
        self.logger.info(f"Fitting {self.model_type} model...")
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=test_size, random_state=random_state
        )
        
        # Scale features
        X_train_scaled = self.scaler.fit_transform(X_train)
        X_test_scaled = self.scaler.transform(X_test)
        
        # Fit model
        self.model.fit(X_train_scaled, y_train)
        
        # Make predictions
        y_train_pred = self.model.predict(X_train_scaled)
        y_test_pred = self.model.predict(X_test_scaled)
        
        # Calculate metrics
        train_mse = mean_squared_error(y_train, y_train_pred)
        test_mse = mean_squared_error(y_test, y_test_pred)
        train_r2 = r2_score(y_train, y_train_pred)
        test_r2 = r2_score(y_test, y_test_pred)
        train_mae = mean_absolute_error(y_train, y_train_pred)
        test_mae = mean_absolute_error(y_test, y_test_pred)
        
        # Cross-validation score
        cv_scores = cross_val_score(self.model, X_train_scaled, y_train, cv=5)
        
        self.is_fitted = True
        
        results = {
            'train_mse': train_mse,
            'test_mse': test_mse,
            'train_r2': train_r2,
            'test_r2': test_r2,
            'train_mae': train_mae,
            'test_mae': test_mae,
            'cv_mean': cv_scores.mean(),
            'cv_std': cv_scores.std(),
            'y_train': y_train,
            'y_test': y_test,
            'y_train_pred': y_train_pred,
            'y_test_pred': y_test_pred
        }
        
        self.logger.info(f"Training completed. Test R²: {test_r2:.4f}")
        
        return results
    
    def predict(self, X: np.ndarray) -> np.ndarray:
        """
        Make predictions.
        
        Parameters:
        -----------
        X : np.ndarray
            Feature matrix
        
        Returns:
        --------
        np.ndarray : Predictions
        """
        if not self.is_fitted:
            raise ValueError("Model must be fitted before prediction")
        
        X_scaled = self.scaler.transform(X)
        return self.model.predict(X_scaled)
    
    def get_feature_importance(self) -> Dict[str, float]:
        """
        Get feature importance scores.
        
        Returns:
        --------
        dict : Feature importance scores
        """
        if not self.is_fitted:
            raise ValueError("Model must be fitted before getting feature importance")
        
        if hasattr(self.model, 'feature_importances_'):
            importance = self.model.feature_importances_
        elif hasattr(self.model, 'coef_'):
            importance = np.abs(self.model.coef_)
        else:
            return {}
        
        return {name: score for name, score in zip(self.feature_names, importance)}
    
    def save_model(self, filepath: str):
        """Save the trained model."""
        model_data = {
            'model': self.model,
            'scaler': self.scaler,
            'model_type': self.model_type,
            'target_property': self.target_property,
            'feature_names': self.feature_names,
            'is_fitted': self.is_fitted
        }
        joblib.dump(model_data, filepath)
        self.logger.info(f"Model saved to {filepath}")
    
    def load_model(self, filepath: str):
        """Load a trained model."""
        model_data = joblib.load(filepath)
        self.model = model_data['model']
        self.scaler = model_data['scaler']
        self.model_type = model_data['model_type']
        self.target_property = model_data['target_property']
        self.feature_names = model_data['feature_names']
        self.is_fitted = model_data['is_fitted']
        self.logger.info(f"Model loaded from {filepath}")
    
    def plot_results(self, results: Dict[str, Any], save_path: str = None):
        """
        Plot training results.
        
        Parameters:
        -----------
        results : dict
            Training results
        save_path : str
            Path to save the plot
        """
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Plot 1: Training vs Test predictions
        axes[0, 0].scatter(results['y_train'], results['y_train_pred'], 
                          alpha=0.6, label='Training', color='blue')
        axes[0, 0].scatter(results['y_test'], results['y_test_pred'], 
                          alpha=0.6, label='Test', color='red')
        axes[0, 0].plot([results['y_train'].min(), results['y_train'].max()], 
                       [results['y_train'].min(), results['y_train'].max()], 
                       'k--', alpha=0.8)
        axes[0, 0].set_xlabel('True Values')
        axes[0, 0].set_ylabel('Predictions')
        axes[0, 0].set_title('Predictions vs True Values')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # Plot 2: Residuals
        train_residuals = results['y_train'] - results['y_train_pred']
        test_residuals = results['y_test'] - results['y_test_pred']
        
        axes[0, 1].scatter(results['y_train_pred'], train_residuals, 
                          alpha=0.6, label='Training', color='blue')
        axes[0, 1].scatter(results['y_test_pred'], test_residuals, 
                          alpha=0.6, label='Test', color='red')
        axes[0, 1].axhline(y=0, color='k', linestyle='--', alpha=0.8)
        axes[0, 1].set_xlabel('Predictions')
        axes[0, 1].set_ylabel('Residuals')
        axes[0, 1].set_title('Residual Plot')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        # Plot 3: Feature importance
        if hasattr(self.model, 'feature_importances_'):
            importance = self.get_feature_importance()
            if importance:
                features = list(importance.keys())
                scores = list(importance.values())
                
                # Sort by importance
                sorted_indices = np.argsort(scores)[::-1]
                features = [features[i] for i in sorted_indices]
                scores = [scores[i] for i in sorted_indices]
                
                axes[1, 0].barh(range(len(features)), scores)
                axes[1, 0].set_yticks(range(len(features)))
                axes[1, 0].set_yticklabels(features)
                axes[1, 0].set_xlabel('Importance')
                axes[1, 0].set_title('Feature Importance')
                axes[1, 0].grid(True, alpha=0.3)
        
        # Plot 4: Metrics summary
        metrics_text = f"""
        Model: {self.model_type}
        Target: {self.target_property}
        
        Training Metrics:
        R² = {results['train_r2']:.4f}
        MSE = {results['train_mse']:.4f}
        MAE = {results['train_mae']:.4f}
        
        Test Metrics:
        R² = {results['test_r2']:.4f}
        MSE = {results['test_mse']:.4f}
        MAE = {results['test_mae']:.4f}
        
        Cross-validation:
        Mean = {results['cv_mean']:.4f}
        Std = {results['cv_std']:.4f}
        """
        axes[1, 1].text(0.1, 0.5, metrics_text, transform=axes[1, 1].transAxes,
                       fontsize=10, verticalalignment='center',
                       bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
        axes[1, 1].set_title('Model Performance')
        axes[1, 1].axis('off')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            self.logger.info(f"Plot saved to {save_path}")
        
        plt.show()


class MDPropertyPredictor(PropertyPredictor):
    """
    Specialized property predictor for molecular dynamics data.
    """
    
    def __init__(self, **kwargs):
        """Initialize MD-specific property predictor."""
        super().__init__(**kwargs)
    
    def predict_diffusion_coefficient(self, 
                                    gro_file: str,
                                    xtc_file: str,
                                    residues: Dict[str, List[str]],
                                    start_frame: int = 0,
                                    stop_frame: int = None,
                                    step_frame: int = 1) -> float:
        """
        Predict diffusion coefficient from trajectory.
        
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
        float : Predicted diffusion coefficient
        """
        # Extract features from trajectory
        trajectory_data = self._load_trajectory_data(
            gro_file, xtc_file, residues, start_frame, stop_frame, step_frame
        )
        
        features = self.extract_features(trajectory_data)
        
        # Make prediction
        prediction = self.predict(features)
        
        return np.mean(prediction)
    
    def _load_trajectory_data(self, 
                             gro_file: str,
                             xtc_file: str,
                             residues: Dict[str, List[str]],
                             start_frame: int,
                             stop_frame: int,
                             step_frame: int) -> Dict[str, np.ndarray]:
        """Load trajectory data from files."""
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
            
        except ImportError:
            self.logger.error("MDAnalysis not available")
            raise
        except Exception as e:
            self.logger.error(f"Error loading trajectory: {e}")
            raise


# Example usage
if __name__ == "__main__":
    # Example: Predict diffusion coefficient
    predictor = MDPropertyPredictor(
        model_type='random_forest',
        target_property='diffusion_coefficient',
        n_estimators=100
    )
    
    # This would require training data first
    # For demonstration, we'll create synthetic data
    np.random.seed(42)
    n_samples = 1000
    n_features = 10
    
    X = np.random.randn(n_samples, n_features)
    y = np.random.randn(n_samples)  # Synthetic target values
    
    # Train model
    results = predictor.fit(X, y)
    
    # Plot results
    predictor.plot_results(results, save_path="property_prediction_results.png")
    
    print(f"Test R²: {results['test_r2']:.4f}")
    print(f"Test MAE: {results['test_mae']:.4f}")
