#!/usr/bin/env python3
"""
Test script for LNB-MDT Machine Learning Module
===============================================

This script tests the basic functionality of the ML module.
"""

import sys
import os
import numpy as np
import pandas as pd
from typing import Dict, List

def test_imports():
    """Test if all ML modules can be imported."""
    print("🔍 Testing imports...")
    
    try:
        from ml import ParameterOptimizer, AnomalyDetector, PropertyPredictor
        from ml import AnalysisParameterOptimizer, MDAnomalyDetector, MDPropertyPredictor
        print("✅ All ML modules imported successfully!")
        return True
    except ImportError as e:
        print(f"❌ Import error: {e}")
        return False

def test_parameter_optimizer():
    """Test parameter optimizer functionality."""
    print("\n🔧 Testing Parameter Optimizer...")
    
    try:
        from ml import ParameterOptimizer
        
        # Define simple parameter bounds
        bounds = {
            'param1': (0, 10),
            'param2': (0, 5)
        }
        
        # Define simple objective function
        def objective_function(params):
            return params['param1']**2 + params['param2']**2
        
        # Create optimizer
        optimizer = ParameterOptimizer(
            parameter_bounds=bounds,
            objective_function=objective_function,
            n_initial_points=5,
            n_iterations=10
        )
        
        # Run optimization
        results = optimizer.optimize()
        
        print(f"✅ Optimization completed!")
        print(f"   Best parameters: {results['best_parameters']}")
        print(f"   Best score: {results['best_score']:.4f}")
        
        return True
        
    except Exception as e:
        print(f"❌ Parameter optimizer test failed: {e}")
        return False

def test_anomaly_detector():
    """Test anomaly detector functionality."""
    print("\n🔍 Testing Anomaly Detector...")
    
    try:
        from ml import AnomalyDetector
        
        # Create synthetic data
        np.random.seed(42)
        n_samples = 100
        n_features = 5
        
        # Generate normal data
        normal_data = np.random.randn(n_samples, n_features)
        
        # Add some anomalies
        anomaly_data = np.random.randn(10, n_features) * 5  # Outliers
        data = np.vstack([normal_data, anomaly_data])
        
        # Create detector
        detector = AnomalyDetector(
            method='isolation_forest',
            contamination=0.1
        )
        
        # Fit and predict
        detector.fit(data)
        predictions = detector.predict(data)
        probabilities = detector.predict_proba(data)
        
        # Get results
        results = detector.detect_anomalies(data)
        
        print(f"✅ Anomaly detection completed!")
        print(f"   Total samples: {len(data)}")
        print(f"   Anomalies detected: {results['n_anomalies']}")
        print(f"   Anomaly ratio: {results['anomaly_ratio']:.2%}")
        
        return True
        
    except Exception as e:
        print(f"❌ Anomaly detector test failed: {e}")
        return False

def test_property_predictor():
    """Test property predictor functionality."""
    print("\n📊 Testing Property Predictor...")
    
    try:
        from ml import PropertyPredictor
        
        # Create synthetic data
        np.random.seed(42)
        n_samples = 200
        n_features = 8
        
        # Generate features and targets
        X = np.random.randn(n_samples, n_features)
        y = np.random.randn(n_samples) * 0.1 + 1.0  # Synthetic property values
        
        # Create predictor
        predictor = PropertyPredictor(
            model_type='random_forest',
            target_property='test_property',
            n_estimators=50,
            max_depth=5
        )
        
        # Train model
        results = predictor.fit(X, y, test_size=0.3)
        
        # Make predictions
        predictions = predictor.predict(X[:10])
        
        # Get feature importance
        importance = predictor.get_feature_importance()
        
        print(f"✅ Property prediction completed!")
        print(f"   Training R²: {results['train_r2']:.4f}")
        print(f"   Test R²: {results['test_r2']:.4f}")
        print(f"   Cross-validation mean: {results['cv_mean']:.4f}")
        print(f"   Number of features: {len(importance)}")
        
        return True
        
    except Exception as e:
        print(f"❌ Property predictor test failed: {e}")
        return False

def test_analysis_parameter_optimizer():
    """Test analysis-specific parameter optimizer."""
    print("\n🎯 Testing Analysis Parameter Optimizer...")
    
    try:
        from ml import AnalysisParameterOptimizer
        
        # Create optimizer for area analysis
        optimizer = AnalysisParameterOptimizer('area')
        
        # Check if bounds are set correctly
        bounds = optimizer.parameter_bounds
        print(f"✅ Analysis parameter optimizer created!")
        print(f"   Analysis type: {optimizer.analysis_type}")
        print(f"   Parameter bounds: {bounds}")
        
        return True
        
    except Exception as e:
        print(f"❌ Analysis parameter optimizer test failed: {e}")
        return False

def test_md_anomaly_detector():
    """Test MD-specific anomaly detector."""
    print("\n🧬 Testing MD Anomaly Detector...")
    
    try:
        from ml import MDAnomalyDetector
        
        # Create detector
        detector = MDAnomalyDetector(
            method='isolation_forest',
            contamination=0.1
        )
        
        # Create synthetic trajectory data
        n_frames = 50
        n_atoms = 20
        
        trajectory_data = {
            'positions': np.random.randn(n_frames, n_atoms, 3),
            'velocities': np.random.randn(n_frames, n_atoms, 3)
        }
        
        # Extract features
        features = detector.extract_features(trajectory_data)
        
        # Fit and detect anomalies
        detector.fit(features)
        results = detector.detect_anomalies(features)
        
        print(f"✅ MD anomaly detection completed!")
        print(f"   Frames analyzed: {len(features)}")
        print(f"   Features extracted: {features.shape[1]}")
        print(f"   Anomalies detected: {results['n_anomalies']}")
        
        return True
        
    except Exception as e:
        print(f"❌ MD anomaly detector test failed: {e}")
        return False

def test_md_property_predictor():
    """Test MD-specific property predictor."""
    print("\n🧬 Testing MD Property Predictor...")
    
    try:
        from ml import MDPropertyPredictor
        
        # Create predictor
        predictor = MDPropertyPredictor(
            model_type='random_forest',
            target_property='diffusion_coefficient',
            n_estimators=50
        )
        
        # Create synthetic trajectory data
        n_frames = 30
        n_atoms = 15
        
        trajectory_data = {
            'positions': np.random.randn(n_frames, n_atoms, 3),
            'velocities': np.random.randn(n_frames, n_atoms, 3)
        }
        
        # Extract features
        features = predictor.extract_features(trajectory_data)
        
        # Create synthetic targets
        targets = np.random.randn(len(features)) * 0.1 + 1.0
        
        # Train model
        results = predictor.fit(features, targets, test_size=0.3)
        
        print(f"✅ MD property prediction completed!")
        print(f"   Features extracted: {features.shape[1]}")
        print(f"   Training R²: {results['train_r2']:.4f}")
        print(f"   Test R²: {results['test_r2']:.4f}")
        
        return True
        
    except Exception as e:
        print(f"❌ MD property predictor test failed: {e}")
        return False

def main():
    """Run all tests."""
    print("🚀 LNB-MDT Machine Learning Module Test Suite")
    print("=" * 50)
    
    tests = [
        ("Import Test", test_imports),
        ("Parameter Optimizer", test_parameter_optimizer),
        ("Anomaly Detector", test_anomaly_detector),
        ("Property Predictor", test_property_predictor),
        ("Analysis Parameter Optimizer", test_analysis_parameter_optimizer),
        ("MD Anomaly Detector", test_md_anomaly_detector),
        ("MD Property Predictor", test_md_property_predictor)
    ]
    
    passed = 0
    total = len(tests)
    
    for test_name, test_func in tests:
        try:
            if test_func():
                passed += 1
        except Exception as e:
            print(f"❌ {test_name} failed with exception: {e}")
    
    print("\n" + "=" * 50)
    print(f"📊 Test Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed! ML module is working correctly.")
        return 0
    else:
        print("⚠️  Some tests failed. Please check the errors above.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
