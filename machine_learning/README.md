# Machine Learning Module for LNB-MDT

This module provides machine learning capabilities for lipid nanobubble analysis.

## 📁 Module Structure

```
machine_learning/
├── __init__.py                 # Module initialization
├── parameter_optimizer.py      # Parameter optimization
├── anomaly_detector.py         # Anomaly detection
├── predictor.py               # Property prediction
├── k_value_optimizer.py       # K-value optimization
├── pattern_recognizer.py      # Pattern recognition
├── data_processor.py          # Data processing
├── README.md                  # This documentation
├── ML_FEATURES_SUMMARY.md     # Feature summary
├── ml_demo.py                 # Interactive demos
├── k_value_demo.py            # K-value optimization demos
├── k_value_example.py         # Simple k-value example
├── test_ml_module.py          # Unit tests
└── test_height_k_value.py     # Height analysis k-value test
```

## 📋 Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Modules](#modules)
- [Examples](#examples)
- [API Reference](#api-reference)

## ✨ Features

### 🔧 Automatic Parameter Optimization
- **Bayesian Optimization**: Uses Gaussian Process regression for efficient parameter search
- **Multi-objective Optimization**: Optimize multiple objectives simultaneously
- **Analysis-specific Bounds**: Predefined parameter ranges for different analysis types
- **Parallel Evaluation**: Support for parallel objective function evaluation

### 🔍 Anomaly Detection
- **Multiple Algorithms**: Isolation Forest, Local Outlier Factor, Elliptic Envelope
- **Feature Extraction**: Automatic extraction of structural, dynamical, and thermodynamic features
- **Visualization**: Comprehensive plotting of anomaly detection results
- **Real-time Analysis**: Detect anomalies in molecular dynamics trajectories

### 📊 Property Prediction
- **Multiple Models**: Random Forest, Gradient Boosting, Neural Networks, SVM, Linear Models
- **Feature Engineering**: Advanced feature extraction from trajectory data
- **Model Evaluation**: Comprehensive metrics and cross-validation
- **Model Persistence**: Save and load trained models

## 🚀 Installation

### Dependencies

The ML module requires additional dependencies beyond the base LNB-MDT requirements:

```bash
pip install scikit-learn scipy matplotlib seaborn joblib
```

### Verify Installation

```python
from ml import ParameterOptimizer, AnomalyDetector, PropertyPredictor
print("ML module installed successfully!")
```

## 🎯 Quick Start

### 1. Parameter Optimization

```python
from ml import AnalysisParameterOptimizer

# Create optimizer for PCA analysis
optimizer = AnalysisParameterOptimizer('pca')

# Create objective function
def objective_function(params):
    # Your analysis function here
    return computation_time + accuracy_penalty

# Run optimization
results = optimizer.optimize()
print(f"Best parameters: {results['best_parameters']}")
```

### 2. Anomaly Detection

```python
from ml import MDAnomalyDetector

# Create anomaly detector
detector = MDAnomalyDetector(method='isolation_forest', contamination=0.1)

# Analyze trajectory
results = detector.analyze_trajectory(
    gro_file="cases/lnb.gro",
    xtc_file="cases/md.xtc",
    residues={'DPPC': ['PO4'], 'CHOL': ['ROH']},
    start_frame=0,
    stop_frame=100
)

# Plot results
detector.plot_anomalies(results, save_path="anomalies.png")
```

### 3. Property Prediction

```python
from ml import MDPropertyPredictor

# Create predictor
predictor = MDPropertyPredictor(
    model_type='random_forest',
    target_property='diffusion_coefficient',
    n_estimators=100
)

# Train model (requires training data)
results = predictor.fit(X_train, y_train)

# Make predictions
predictions = predictor.predict(X_test)

# Plot results
predictor.plot_results(results, save_path="predictions.png")
```

## 📦 Modules

### ParameterOptimizer

Automatic parameter optimization using Bayesian optimization.

**Key Features:**
- Gaussian Process-based optimization
- Support for custom objective functions
- Analysis-specific parameter bounds
- Parallel evaluation support

**Example:**
```python
from ml import ParameterOptimizer

# Define parameter bounds
bounds = {
    'k_value': (5, 50),
    'cutoff': (5.0, 15.0),
    'n_components': (2, 10)
}

# Create optimizer
optimizer = ParameterOptimizer(bounds, objective_function)

# Run optimization
results = optimizer.optimize()
```

### AnomalyDetector

Detect anomalies in molecular dynamics trajectories.

**Key Features:**
- Multiple detection algorithms
- Automatic feature extraction
- Comprehensive visualization
- Real-time analysis

**Example:**
```python
from ml import AnomalyDetector

# Create detector
detector = AnomalyDetector(method='isolation_forest', contamination=0.1)

# Fit and predict
detector.fit(features)
predictions = detector.predict(new_features)
probabilities = detector.predict_proba(new_features)
```

### PropertyPredictor

Predict molecular properties from trajectory data.

**Key Features:**
- Multiple ML algorithms
- Advanced feature engineering
- Model evaluation and visualization
- Model persistence

**Example:**
```python
from ml import PropertyPredictor

# Create predictor
predictor = PropertyPredictor(
    model_type='random_forest',
    target_property='diffusion_coefficient'
)

# Train model
results = predictor.fit(X, y)

# Get feature importance
importance = predictor.get_feature_importance()

# Save model
predictor.save_model('model.pkl')
```

## 📊 Examples

### Complete Parameter Optimization Workflow

```python
from ml import AnalysisParameterOptimizer
import time

# Create optimizer for area analysis
optimizer = AnalysisParameterOptimizer('area')

# Define objective function
def objective_function(params):
    try:
        # Run area analysis with given parameters
        from analysis.area import Area
        
        analyzer = Area(
            gro_file="cases/lnb.gro",
            xtc_file="cases/md.xtc",
            residues={'DPPC': ['PO4']},
            **params
        )
        
        start_time = time.time()
        results = analyzer.run()
        computation_time = time.time() - start_time
        
        # Calculate objective (lower is better)
        # Consider both computation time and result quality
        objective = computation_time + len(results) * 0.001
        
        return objective
        
    except Exception as e:
        print(f"Error: {e}")
        return float('inf')  # High penalty for failed evaluations

# Run optimization
results = optimizer.optimize()

print("Optimization completed!")
print(f"Best parameters: {results['best_parameters']}")
print(f"Best score: {results['best_score']}")

# Save optimized model
optimizer.save_model('optimized_area_params.pkl')
```

### Anomaly Detection in Trajectory

```python
from ml import MDAnomalyDetector
import matplotlib.pyplot as plt

# Create detector
detector = MDAnomalyDetector(
    method='isolation_forest',
    contamination=0.1
)

# Analyze trajectory
results = detector.analyze_trajectory(
    gro_file="cases/lnb.gro",
    xtc_file="cases/md.xtc",
    residues={'DPPC': ['PO4'], 'CHOL': ['ROH']},
    start_frame=0,
    stop_frame=1000,
    step_frame=5
)

# Print results
print(f"Total frames analyzed: {len(results['predictions'])}")
print(f"Anomalies detected: {results['n_anomalies']}")
print(f"Anomaly ratio: {results['anomaly_ratio']:.2%}")

# Plot results
detector.plot_anomalies(results, save_path="anomaly_analysis.png")

# Analyze specific anomalies
anomaly_frames = results['anomaly_indices']
print(f"Anomaly frames: {anomaly_frames}")
```

### Property Prediction Pipeline

```python
from ml import MDPropertyPredictor
import numpy as np
import pandas as pd

# Create predictor
predictor = MDPropertyPredictor(
    model_type='random_forest',
    target_property='diffusion_coefficient',
    n_estimators=100,
    max_depth=10
)

# Generate synthetic training data (replace with real data)
np.random.seed(42)
n_samples = 1000
n_features = 15

# Create synthetic features and targets
X = np.random.randn(n_samples, n_features)
y = np.random.randn(n_samples) * 0.1 + 1.0  # Synthetic diffusion coefficients

# Train model
results = predictor.fit(X, y, test_size=0.2)

# Print performance metrics
print(f"Training R²: {results['train_r2']:.4f}")
print(f"Test R²: {results['test_r2']:.4f}")
print(f"Cross-validation mean: {results['cv_mean']:.4f}")

# Get feature importance
importance = predictor.get_feature_importance()
print("Top 5 most important features:")
for feature, score in sorted(importance.items(), key=lambda x: x[1], reverse=True)[:5]:
    print(f"  {feature}: {score:.4f}")

# Plot results
predictor.plot_results(results, save_path="property_prediction.png")

# Save model
predictor.save_model('diffusion_predictor.pkl')

# Make predictions on new data
new_features = np.random.randn(10, n_features)
predictions = predictor.predict(new_features)
print(f"Predictions: {predictions}")
```

## 🔧 API Reference

### ParameterOptimizer

#### `__init__(parameter_bounds, objective_function, **kwargs)`
Initialize the parameter optimizer.

**Parameters:**
- `parameter_bounds`: Dict mapping parameter names to (min, max) bounds
- `objective_function`: Function that takes parameters and returns objective value
- `n_initial_points`: Number of initial random points (default: 10)
- `n_iterations`: Number of optimization iterations (default: 50)
- `random_state`: Random seed (default: 42)

#### `optimize()`
Run the optimization process.

**Returns:**
- Dict containing best parameters, best score, and optimization history

#### `save_model(filepath)`
Save the trained model to file.

#### `load_model(filepath)`
Load a trained model from file.

### AnomalyDetector

#### `__init__(method='isolation_forest', **kwargs)`
Initialize the anomaly detector.

**Parameters:**
- `method`: Detection method ('isolation_forest', 'lof', 'elliptic_envelope')
- `contamination`: Expected fraction of anomalies (default: 0.1)

#### `fit(data)`
Fit the anomaly detection model.

#### `predict(data)`
Predict anomalies in the data.

**Returns:**
- Array of predictions (-1 for anomalies, 1 for normal)

#### `predict_proba(data)`
Predict anomaly probabilities.

**Returns:**
- Array of anomaly probabilities

#### `detect_anomalies(data, threshold=0.5)`
Detect anomalies with comprehensive results.

**Returns:**
- Dict containing predictions, probabilities, and statistics

### PropertyPredictor

#### `__init__(model_type='random_forest', target_property='diffusion_coefficient', **kwargs)`
Initialize the property predictor.

**Parameters:**
- `model_type`: Type of model ('random_forest', 'gradient_boosting', 'linear', 'svr', 'neural_network')
- `target_property`: Property to predict
- Additional model-specific parameters

#### `fit(X, y, test_size=0.2, random_state=42)`
Fit the model to the data.

**Returns:**
- Dict containing training results and metrics

#### `predict(X)`
Make predictions on new data.

**Returns:**
- Array of predictions

#### `get_feature_importance()`
Get feature importance scores.

**Returns:**
- Dict mapping feature names to importance scores

#### `plot_results(results, save_path=None)`
Plot training results and model performance.

#### `save_model(filepath)`
Save the trained model to file.

#### `load_model(filepath)`
Load a trained model from file.

## 🎯 Best Practices

### Parameter Optimization
1. **Define Clear Objectives**: Ensure your objective function captures both accuracy and efficiency
2. **Set Appropriate Bounds**: Use domain knowledge to set reasonable parameter ranges
3. **Monitor Progress**: Use logging to track optimization progress
4. **Validate Results**: Always validate optimized parameters on test data

### Anomaly Detection
1. **Choose Appropriate Method**: Isolation Forest for general anomalies, LOF for local anomalies
2. **Tune Contamination**: Set contamination based on expected anomaly ratio
3. **Feature Selection**: Use relevant features for your specific analysis
4. **Interpret Results**: Analyze detected anomalies in the context of your system

### Property Prediction
1. **Feature Engineering**: Extract meaningful features from trajectory data
2. **Model Selection**: Try multiple models and compare performance
3. **Cross-validation**: Use cross-validation to assess model generalization
4. **Feature Importance**: Analyze feature importance to understand predictions

## 🚨 Troubleshooting

### Common Issues

1. **Import Errors**
   ```bash
   pip install scikit-learn scipy matplotlib seaborn joblib
   ```

2. **Memory Issues**
   - Reduce batch size for large datasets
   - Use feature selection to reduce dimensionality
   - Process data in chunks

3. **Poor Performance**
   - Check feature scaling
   - Try different model types
   - Adjust hyperparameters
   - Increase training data

4. **Convergence Issues**
   - Check parameter bounds
   - Adjust optimization parameters
   - Verify objective function

## 📚 References

- [Scikit-learn Documentation](https://scikit-learn.org/)
- [Bayesian Optimization](https://arxiv.org/abs/1807.02811)
- [Anomaly Detection in Time Series](https://arxiv.org/abs/1901.03407)
- [Molecular Property Prediction](https://pubs.acs.org/doi/10.1021/acs.jcim.9b00237)

---

**🤖 LNB-MDT Machine Learning Module** - Making molecular dynamics analysis smarter and more efficient!
