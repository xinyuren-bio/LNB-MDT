#!/usr/bin/env python3
"""
K-Value Optimization Demo for LNB-MDT
====================================

This script demonstrates how to automatically find the optimal k-value
for different molecular dynamics analysis types.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List
import time

def demo_k_value_optimization():
    """Demonstrate k-value optimization for different analysis types."""
    print("🎯 K-Value Optimization Demo")
    print("=" * 50)
    
    try:
        from .k_value_optimizer import KValueOptimizer
        
        # Test different analysis types
        analysis_types = ['area', 'curvature', 'height', 'cluster']
        
        results_summary = {}
        
        for analysis_type in analysis_types:
            print(f"\n🔧 Optimizing k-value for {analysis_type} analysis...")
            
            # Create optimizer
            optimizer = KValueOptimizer(
                analysis_type=analysis_type,
                n_trials=10,  # Reduced for demo
                k_range=(5, 40)
            )
            
            # Run optimization
            results = optimizer.optimize(
                gro_file="cases/lnb.gro",
                xtc_file="cases/md.xtc",
                residues={'DPPC': ['PO4'], 'CHOL': ['ROH']}
            )
            
            # Store results
            results_summary[analysis_type] = results
            
            print(f"✅ {analysis_type} optimization completed!")
            print(f"   Optimal k-value: {results['best_k']}")
            print(f"   Best score: {results['best_score']:.4f}")
            print(f"   K-values tested: {len(results['k_values_tested'])}")
            
            # Plot results
            optimizer.plot_optimization_results(
                results, 
                save_path=f"k_value_optimization_{analysis_type}.png"
            )
            
            # Get recommendations
            recommendations = optimizer.get_recommendations(results)
            print(f"   Recommendations: {recommendations['reasoning']}")
        
        # Compare results across analysis types
        print("\n📊 Comparison Across Analysis Types:")
        print("-" * 40)
        for analysis_type, results in results_summary.items():
            print(f"{analysis_type:12}: k={results['best_k']:2d} (score: {results['best_score']:.4f})")
        
        return True
        
    except Exception as e:
        print(f"❌ K-value optimization demo failed: {e}")
        return False

def demo_smart_k_value_selector():
    """Demonstrate smart k-value selector that learns from history."""
    print("\n🧠 Smart K-Value Selector Demo")
    print("=" * 50)
    
    try:
        from .k_value_optimizer import SmartKValueSelector, KValueOptimizer
        
        # Create smart selector
        smart_selector = SmartKValueSelector()
        
        # Simulate optimization history for different systems
        print("Training smart selector with optimization history...")
        
        # System 1: Small system
        system1_features = {
            'n_atoms': 500,
            'density': 0.05,
            'box_length': 8.0,
            'n_residues': 50,
            'avg_residue_size': 10.0
        }
        
        # System 2: Medium system
        system2_features = {
            'n_atoms': 1000,
            'density': 0.1,
            'box_length': 10.0,
            'n_residues': 100,
            'avg_residue_size': 10.0
        }
        
        # System 3: Large system
        system3_features = {
            'n_atoms': 2000,
            'density': 0.15,
            'box_length': 12.0,
            'n_residues': 200,
            'avg_residue_size': 10.0
        }
        
        # Add optimization results to history
        optimization_results = [
            {'system_features': system1_features, 'best_k': 12},
            {'system_features': system2_features, 'best_k': 20},
            {'system_features': system3_features, 'best_k': 28},
            {'system_features': {'n_atoms': 800, 'density': 0.08, 'box_length': 9.0, 'n_residues': 80, 'avg_residue_size': 10.0}, 'best_k': 16},
            {'system_features': {'n_atoms': 1500, 'density': 0.12, 'box_length': 11.0, 'n_residues': 150, 'avg_residue_size': 10.0}, 'best_k': 24}
        ]
        
        for result in optimization_results:
            smart_selector.add_optimization_result(result)
        
        # Train the model
        smart_selector.train_model()
        
        # Test prediction on new systems
        print("\nTesting predictions on new systems:")
        
        test_systems = [
            {'name': 'Small dense system', 'features': {'n_atoms': 600, 'density': 0.12, 'box_length': 7.0, 'n_residues': 60, 'avg_residue_size': 10.0}},
            {'name': 'Large sparse system', 'features': {'n_atoms': 1800, 'density': 0.06, 'box_length': 15.0, 'n_residues': 180, 'avg_residue_size': 10.0}},
            {'name': 'Medium balanced system', 'features': {'n_atoms': 1200, 'density': 0.1, 'box_length': 10.5, 'n_residues': 120, 'avg_residue_size': 10.0}}
        ]
        
        for test_system in test_systems:
            predicted_k = smart_selector.predict_optimal_k(test_system['features'])
            print(f"   {test_system['name']}: predicted k = {predicted_k}")
        
        return True
        
    except Exception as e:
        print(f"❌ Smart k-value selector demo failed: {e}")
        return False

def demo_k_value_impact():
    """Demonstrate the impact of different k-values on analysis results."""
    print("\n📈 K-Value Impact Demo")
    print("=" * 50)
    
    try:
        from .k_value_optimizer import KValueOptimizer
        
        # Test different k-values for area analysis
        k_values = [10, 15, 20, 25, 30, 35]
        computation_times = []
        result_qualities = []
        
        print("Testing different k-values for area analysis...")
        
        for k_value in k_values:
            print(f"   Testing k = {k_value}...")
            
            # Create optimizer to test specific k-value
            optimizer = KValueOptimizer('area')
            
            # Test the k-value
            start_time = time.time()
            score = optimizer._evaluate_k_value(
                k_value, 
                "cases/lnb.gro", 
                "cases/md.xtc",
                {'DPPC': ['PO4'], 'CHOL': ['ROH']}
            )
            computation_time = time.time() - start_time
            
            computation_times.append(computation_time)
            result_qualities.append(1.0 / (1.0 + score))  # Convert score to quality
            
            print(f"     Computation time: {computation_time:.3f}s")
            print(f"     Quality score: {result_qualities[-1]:.4f}")
        
        # Plot impact
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Plot computation time
        ax1.plot(k_values, computation_times, 'bo-', linewidth=2, markersize=8)
        ax1.set_xlabel('K-Value')
        ax1.set_ylabel('Computation Time (s)')
        ax1.set_title('Computation Time vs K-Value')
        ax1.grid(True, alpha=0.3)
        
        # Plot result quality
        ax2.plot(k_values, result_qualities, 'ro-', linewidth=2, markersize=8)
        ax2.set_xlabel('K-Value')
        ax2.set_ylabel('Result Quality')
        ax2.set_title('Result Quality vs K-Value')
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('k_value_impact.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # Find optimal k-value
        optimal_idx = np.argmax(result_qualities)
        optimal_k = k_values[optimal_idx]
        
        print(f"\n📊 Impact Analysis Results:")
        print(f"   Optimal k-value: {optimal_k}")
        print(f"   Best quality: {result_qualities[optimal_idx]:.4f}")
        print(f"   Computation time: {computation_times[optimal_idx]:.3f}s")
        
        return True
        
    except Exception as e:
        print(f"❌ K-value impact demo failed: {e}")
        return False

def demo_practical_workflow():
    """Demonstrate practical workflow for k-value optimization."""
    print("\n🔄 Practical Workflow Demo")
    print("=" * 50)
    
    try:
        from .k_value_optimizer import KValueOptimizer
        
        print("Step 1: Analyze your system properties...")
        optimizer = KValueOptimizer('area')
        features = optimizer._extract_system_features("cases/lnb.gro", "cases/md.xtc")
        
        print(f"   System features:")
        for key, value in features.items():
            print(f"     {key}: {value}")
        
        print("\nStep 2: Get initial k-value estimate...")
        initial_k = optimizer._calculate_optimal_k_heuristic(features)
        print(f"   Initial k-value estimate: {initial_k}")
        
        print("\nStep 3: Run optimization...")
        results = optimizer.optimize(
            gro_file="cases/lnb.gro",
            xtc_file="cases/md.xtc",
            residues={'DPPC': ['PO4'], 'CHOL': ['ROH']}
        )
        
        print(f"   Optimal k-value: {results['best_k']}")
        print(f"   Optimization completed in {len(results['k_values_tested'])} trials")
        
        print("\nStep 4: Get recommendations...")
        recommendations = optimizer.get_recommendations(results)
        print(f"   Confidence: {recommendations['confidence']}")
        print(f"   Reasoning: {recommendations['reasoning']}")
        print(f"   Alternative k-values: {recommendations['alternative_k_values']}")
        
        print("\nStep 5: Use optimal k-value in your analysis...")
        print(f"   python analysis/area.py --k-value {results['best_k']} --gro-file cases/lnb.gro --xtc-file cases/md.xtc")
        
        return True
        
    except Exception as e:
        print(f"❌ Practical workflow demo failed: {e}")
        return False

def main():
    """Run all k-value optimization demos."""
    print("🚀 LNB-MDT K-Value Optimization Demo")
    print("=" * 60)
    print("This demo shows how to automatically find optimal k-values")
    print("for different molecular dynamics analysis types.")
    print()
    
    demos = [
        ("K-Value Optimization", demo_k_value_optimization),
        ("Smart K-Value Selector", demo_smart_k_value_selector),
        ("K-Value Impact Analysis", demo_k_value_impact),
        ("Practical Workflow", demo_practical_workflow)
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
    
    print("=" * 60)
    print(f"📊 Demo Results: {passed}/{total} demos completed successfully")
    
    if passed == total:
        print("🎉 All demos completed! K-value optimization is working correctly.")
        print("\nKey benefits:")
        print("1. Automatic k-value selection saves hours of manual tuning")
        print("2. Optimal k-values improve analysis accuracy")
        print("3. System-specific recommendations based on properties")
        print("4. Learning from history for future predictions")
    else:
        print("⚠️  Some demos failed. Please check the errors above.")
    
    return 0 if passed == total else 1

if __name__ == "__main__":
    import sys
    sys.exit(main())
