"""
Parameter Optimizer for LNB-MDT
===============================

This module provides automatic parameter optimization using machine learning techniques.
"""

import numpy as np
import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel
from sklearn.model_selection import cross_val_score
from scipy.optimize import minimize
import joblib
from typing import Dict, List, Tuple, Any, Callable
import logging
import time
from scipy import stats as norm

class ParameterOptimizer:
    """
    Automatic parameter optimizer using Bayesian optimization.
    
    This class implements Bayesian optimization to automatically find optimal
    parameters for molecular dynamics analysis.
    """
    
    def __init__(self, 
                 parameter_bounds: Dict[str, Tuple[float, float]],
                 objective_function: Callable,
                 n_initial_points: int = 10,
                 n_iterations: int = 50,
                 random_state: int = 42):
        """
        Initialize the parameter optimizer.
        
        Parameters:
        -----------
        parameter_bounds : dict
            Dictionary mapping parameter names to (min, max) bounds
        objective_function : callable
            Function that takes parameters and returns objective value
        n_initial_points : int
            Number of initial random points to sample
        n_iterations : int
            Number of optimization iterations
        random_state : int
            Random seed for reproducibility
        """
        self.parameter_bounds = parameter_bounds
        self.objective_function = objective_function
        self.n_initial_points = n_initial_points
        self.n_iterations = n_iterations
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
        self.X_history = []
        self.y_history = []
        self.best_params = None
        self.best_score = float('inf')
        
        # Setup logging
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)
    
    def _sample_random_points(self, n_points: int) -> np.ndarray:
        """Sample random points within parameter bounds."""
        points = []
        for _ in range(n_points):
            point = []
            for param_name, (min_val, max_val) in self.parameter_bounds.items():
                point.append(np.random.uniform(min_val, max_val))
            points.append(point)
        return np.array(points)
    
    def _acquisition_function(self, X: np.ndarray) -> np.ndarray:
        """Expected Improvement acquisition function."""
        if len(self.X_history) == 0:
            return np.zeros(len(X))
        
        X = X.reshape(-1, len(self.parameter_bounds))
        mean, std = self.gp.predict(X, return_std=True)
        
        # Expected Improvement
        best_f = min(self.y_history)
        improvement = best_f - mean
        z = improvement / (std + 1e-8)
        ei = improvement * norm.cdf(z) + std * norm.pdf(z)
        
        return ei
    
    def optimize(self) -> Dict[str, Any]:
        """
        Run the optimization process.
        
        Returns:
        --------
        dict : Optimization results including best parameters and history
        """
        self.logger.info("Starting parameter optimization...")
        
        # Initial random sampling
        X_initial = self._sample_random_points(self.n_initial_points)
        
        for i, x in enumerate(X_initial):
            self.logger.info(f"Evaluating initial point {i+1}/{self.n_initial_points}")
            y = self.objective_function(self._params_dict(x))
            self.X_history.append(x)
            self.y_history.append(y)
            
            if y < self.best_score:
                self.best_score = y
                self.best_params = self._params_dict(x)
        
        # Bayesian optimization iterations
        for iteration in range(self.n_iterations):
            self.logger.info(f"Optimization iteration {iteration+1}/{self.n_iterations}")
            
            # Fit Gaussian Process
            X_train = np.array(self.X_history)
            y_train = np.array(self.y_history)
            self.gp.fit(X_train, y_train)
            
            # Find next point to evaluate
            next_point = self._find_next_point()
            y_next = self.objective_function(self._params_dict(next_point))
            
            self.X_history.append(next_point)
            self.y_history.append(y_next)
            
            if y_next < self.best_score:
                self.best_score = y_next
                self.best_params = self._params_dict(next_point)
                self.logger.info(f"New best score: {self.best_score}")
        
        return {
            'best_parameters': self.best_params,
            'best_score': self.best_score,
            'optimization_history': {
                'parameters': self.X_history,
                'scores': self.y_history
            }
        }
    
    def _params_dict(self, x: np.ndarray) -> Dict[str, float]:
        """Convert parameter array to dictionary."""
        return {name: val for name, val in zip(self.parameter_bounds.keys(), x)}
    
    def _find_next_point(self) -> np.ndarray:
        """Find the next point to evaluate using acquisition function."""
        # Grid search over parameter space
        n_grid = 1000
        X_grid = self._sample_random_points(n_grid)
        
        # Evaluate acquisition function
        acq_values = self._acquisition_function(X_grid)
        
        # Return point with maximum acquisition value
        return X_grid[np.argmax(acq_values)]
    
    def save_model(self, filepath: str):
        """Save the trained model."""
        model_data = {
            'gp': self.gp,
            'parameter_bounds': self.parameter_bounds,
            'best_params': self.best_params,
            'best_score': self.best_score,
            'X_history': self.X_history,
            'y_history': self.y_history
        }
        joblib.dump(model_data, filepath)
        self.logger.info(f"Model saved to {filepath}")
    
    def load_model(self, filepath: str):
        """Load a trained model."""
        model_data = joblib.load(filepath)
        self.gp = model_data['gp']
        self.parameter_bounds = model_data['parameter_bounds']
        self.best_params = model_data['best_params']
        self.best_score = model_data['best_score']
        self.X_history = model_data['X_history']
        self.y_history = model_data['y_history']
        self.logger.info(f"Model loaded from {filepath}")


class AnalysisParameterOptimizer(ParameterOptimizer):
    """
    Specialized parameter optimizer for LNB-MDT analysis parameters.
    """
    
    def __init__(self, analysis_type: str, **kwargs):
        """
        Initialize with predefined parameter bounds for specific analysis types.
        
        Parameters:
        -----------
        analysis_type : str
            Type of analysis ('pca', 'area', 'curvature', 'cluster', etc.)
        """
        self.analysis_type = analysis_type
        parameter_bounds = self._get_default_bounds(analysis_type)
        super().__init__(parameter_bounds, **kwargs)
    
    def _get_default_bounds(self, analysis_type: str) -> Dict[str, Tuple[float, float]]:
        """Get default parameter bounds for different analysis types."""
        bounds = {
            'pca': {
                'n_components': (2, 10),
                'random_state': (0, 100)
            },
            'area': {
                'k_value': (5, 50),
                'max_normal_angle': (90, 180)
            },
            'curvature': {
                'k_value': (5, 50),
                'smoothing_factor': (0.1, 2.0)
            },
            'cluster': {
                'cutoff': (5.0, 15.0),
                'min_cluster_size': (2, 20)
            },
            'height': {
                'k_value': (5, 50),
                'reference_height': (-10, 10)
            }
        }
        return bounds.get(analysis_type, {})
    
    def create_objective_function(self, 
                                gro_file: str, 
                                xtc_file: str,
                                residues: Dict[str, List[str]]) -> Callable:
        """
        Create an objective function for the specific analysis type.
        
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
        callable : Objective function
        """
        def objective_function(params):
            try:
                # Import analysis modules dynamically
                if self.analysis_type == 'pca':
                    from analysis.pca import PCA
                    analyzer = PCA(gro_file, xtc_file, residues, **params)
                elif self.analysis_type == 'area':
                    from analysis.area import Area
                    analyzer = Area(gro_file, xtc_file, residues, **params)
                elif self.analysis_type == 'curvature':
                    from analysis.curvature import Curvature
                    analyzer = Curvature(gro_file, xtc_file, residues, **params)
                elif self.analysis_type == 'cluster':
                    from analysis.cluster import Cluster
                    analyzer = Cluster(gro_file, xtc_file, residues, **params)
                elif self.analysis_type == 'height':
                    from analysis.height import Height
                    analyzer = Height(gro_file, xtc_file, residues, **params)
                else:
                    raise ValueError(f"Unsupported analysis type: {self.analysis_type}")
                
                # Run analysis and calculate objective (e.g., computation time, accuracy)
                start_time = time.time()
                results = analyzer.run()
                computation_time = time.time() - start_time
                
                # Calculate objective value (lower is better)
                # This is a simple example - you might want to use more sophisticated metrics
                objective = computation_time + len(results) * 0.001  # Penalty for large results
                
                return objective
                
            except Exception as e:
                self.logger.warning(f"Error in objective function: {e}")
                return float('inf')  # Return high penalty for failed evaluations
        
        return objective_function


# Example usage
if __name__ == "__main__":
    # Example: Optimize PCA parameters
    optimizer = AnalysisParameterOptimizer('pca')
    
    # Create objective function
    objective_func = optimizer.create_objective_function(
        gro_file="cases/lnb.gro",
        xtc_file="cases/md.xtc",
        residues={'DPPC': ['PO4'], 'CHOL': ['ROH']}
    )
    
    # Run optimization
    results = optimizer.optimize()
    
    print("Best parameters:", results['best_parameters'])
    print("Best score:", results['best_score'])
