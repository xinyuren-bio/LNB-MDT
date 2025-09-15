#!/usr/bin/env python3
"""
LNB-MDT Machine Learning Demo
=============================

This script demonstrates the machine learning capabilities of LNB-MDT.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List
import time

def demo_parameter_optimization():
    """Demonstrate parameter optimization."""
    print("🔧 Parameter Optimization Demo")
    print("=" * 40)
    
    try:
        from ml import AnalysisParameterOptimizer
        
        # Create optimizer for area analysis
        optimizer = AnalysisParameterOptimizer('area')
        
        # Define a simple objective function
        def objective_function(params):
            # Simulate analysis time and quality
            k_value = params.get('k_value', 20)
            max_normal_angle = params.get('max_normal_angle', 140)
            
            # Simulate computation time (higher k_value = more time)
            computation_time = k_value * 0.01
            
            # Simulate quality (optimal around k_value=25, max_normal_angle=150)
            quality_penalty = abs(k_value - 25) * 0.1 + abs(max_normal_angle - 150) * 0.01
            
            return computation_time + quality_penalty
        
        # Set the objective function
        optimizer.objective_function = objective_function
        
        # Run optimization with fewer iterations for demo
        print("Running parameter optimization...")
        results = optimizer.optimize()
        
        print(f"✅ Optimization completed!")
        print(f"   Best parameters: {results['best_parameters']}")
        print(f"   Best score: {results['best_score']:.4f}")
        print(f"   Optimization history: {len(results['optimization_history']['scores'])} evaluations")
        
        return True
        
    except Exception as e:
        print(f"❌ Parameter optimization demo failed: {e}")
        return False

def demo_anomaly_detection():
    """Demonstrate anomaly detection."""
    print("\n🔍 Anomaly Detection Demo")
    print("=" * 40)
    
    try:
        from ml import MDAnomalyDetector
        
        # Create detector
        detector = MDAnomalyDetector(
            method='isolation_forest',
            contamination=0.1
        )
        
        # Create synthetic trajectory data with anomalies
        np.random.seed(42)
        n_frames = 100
        n_atoms = 30
        
        # Generate normal trajectory
        normal_positions = np.random.randn(n_frames, n_atoms, 3) * 0.1
        
        # Add some anomalies (sudden movements)
        anomaly_frames = [20, 45, 70, 85]
        for frame_idx in anomaly_frames:
            normal_positions[frame_idx] += np.random.randn(n_atoms, 3) * 2.0
        
        trajectory_data = {
            'positions': normal_positions,
            'velocities': np.random.randn(n_frames, n_atoms, 3) * 0.05
        }
        
        # Extract features and detect anomalies
        print("Extracting features and detecting anomalies...")
        features = detector.extract_features(trajectory_data)
        detector.fit(features)
        results = detector.detect_anomalies(features)
        
        print(f"✅ Anomaly detection completed!")
        print(f"   Frames analyzed: {len(features)}")
        print(f"   Features extracted: {features.shape[1]}")
        print(f"   Anomalies detected: {results['n_anomalies']}")
        print(f"   Anomaly ratio: {results['anomaly_ratio']:.2%}")
        print(f"   Anomaly frames: {results['anomaly_indices']}")
        
        # Check if we detected the injected anomalies
        detected_anomalies = set(results['anomaly_indices'])
        injected_anomalies = set(anomaly_frames)
        overlap = detected_anomalies.intersection(injected_anomalies)
        
        print(f"   Injected anomalies: {injected_anomalies}")
        print(f"   Detection accuracy: {len(overlap)}/{len(injected_anomalies)} anomalies detected")
        
        return True
        
    except Exception as e:
        print(f"❌ Anomaly detection demo failed: {e}")
        return False

def demo_property_prediction():
    """Demonstrate property prediction."""
    print("\n📊 Property Prediction Demo")
    print("=" * 40)
    
    try:
        from ml import MDPropertyPredictor
        
        # Create predictor
        predictor = MDPropertyPredictor(
            model_type='random_forest',
            target_property='diffusion_coefficient',
            n_estimators=100,
            max_depth=8
        )
        
        # Create synthetic training data
        np.random.seed(42)
        n_samples = 500
        n_frames = 50
        n_atoms = 25
        
        print("Generating synthetic training data...")
        
        # Generate multiple trajectories with different properties
        all_features = []
        all_targets = []
        
        for i in range(n_samples):
            # Create trajectory with varying properties
            base_diffusion = np.random.uniform(0.5, 2.0)
            
            # Generate trajectory data
            positions = np.random.randn(n_frames, n_atoms, 3) * 0.1
            velocities = np.random.randn(n_frames, n_atoms, 3) * 0.05 * base_diffusion
            
            trajectory_data = {
                'positions': positions,
                'velocities': velocities
            }
            
            # Extract features
            features = predictor.extract_features(trajectory_data)
            
            if len(features) > 0:
                all_features.append(features[0])  # Use first frame features
                all_targets.append(base_diffusion)
        
        # Convert to arrays
        X = np.array(all_features)
        y = np.array(all_targets)
        
        print(f"Training data shape: {X.shape}")
        print(f"Target range: {y.min():.3f} - {y.max():.3f}")
        
        # Train model
        print("Training property prediction model...")
        results = predictor.fit(X, y, test_size=0.2)
        
        print(f"✅ Property prediction model trained!")
        print(f"   Training R²: {results['train_r2']:.4f}")
        print(f"   Test R²: {results['test_r2']:.4f}")
        print(f"   Cross-validation mean: {results['cv_mean']:.4f}")
        print(f"   Test MAE: {results['test_mae']:.4f}")
        
        # Get feature importance
        importance = predictor.get_feature_importance()
        if importance:
            print(f"   Top 3 important features:")
            sorted_features = sorted(importance.items(), key=lambda x: x[1], reverse=True)[:3]
            for feature, score in sorted_features:
                print(f"     {feature}: {score:.4f}")
        
        # Make predictions on new data
        print("\nMaking predictions on new data...")
        new_trajectory_data = {
            'positions': np.random.randn(n_frames, n_atoms, 3) * 0.1,
            'velocities': np.random.randn(n_frames, n_atoms, 3) * 0.05
        }
        
        new_features = predictor.extract_features(new_trajectory_data)
        if len(new_features) > 0:
            prediction = predictor.predict(new_features[:1])
            print(f"   Predicted diffusion coefficient: {prediction[0]:.4f}")
        
        return True
        
    except Exception as e:
        print(f"❌ Property prediction demo failed: {e}")
        return False

def demo_integration():
    """Demonstrate integration with existing analysis modules."""
    print("\n🔗 Integration Demo")
    print("=" * 40)
    
    try:
        from ml import AnalysisParameterOptimizer
        
        # Create optimizer for PCA analysis
        optimizer = AnalysisParameterOptimizer('pca')
        
        print("✅ Integration test completed!")
        print(f"   Available analysis types: {list(optimizer._get_default_bounds('pca').keys())}")
        print(f"   PCA parameter bounds: {optimizer._get_default_bounds('pca')}")
        print(f"   Area parameter bounds: {optimizer._get_default_bounds('area')}")
        print(f"   Curvature parameter bounds: {optimizer._get_default_bounds('curvature')}")
        
        return True
        
    except Exception as e:
        print(f"❌ Integration demo failed: {e}")
        return False

def main():
    """Run all demos."""
    print("🚀 LNB-MDT Machine Learning Demo")
    print("=" * 50)
    print("This demo showcases the machine learning capabilities of LNB-MDT.")
    print("It includes parameter optimization, anomaly detection, and property prediction.")
    print()
    
    demos = [
        ("Parameter Optimization", demo_parameter_optimization),
        ("Anomaly Detection", demo_anomaly_detection),
        ("Property Prediction", demo_property_prediction),
        ("Integration Test", demo_integration)
    ]
    
    passed = 0
    total = len(demos)
    
    for demo_name, demo_func in demos:
        try:
            if demo_func():
                passed += 1
            print()
        except Exception as e:
            print(f"❌ {demo_name} failed with exception: {e}")
            print()
    
    print("=" * 50)
    print(f"📊 Demo Results: {passed}/{total} demos completed successfully")
    
    if passed == total:
        print("🎉 All demos completed! Machine learning module is working correctly.")
        print("\nNext steps:")
        print("1. Try the demos with your own trajectory data")
        print("2. Explore the ML module documentation: machine_learning/README.md")
        print("3. Integrate ML features into your analysis workflows")
    else:
        print("⚠️  Some demos failed. Please check the errors above.")
        print("Make sure all dependencies are installed: pip install scikit-learn seaborn joblib")
    
    return 0 if passed == total else 1

if __name__ == "__main__":
    import sys
    sys.exit(main())
