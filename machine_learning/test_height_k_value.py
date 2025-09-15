#!/usr/bin/env python3
"""
Test K-Value Optimization for Height Analysis
============================================

This script specifically tests k-value optimization for height analysis.
"""

import time
import numpy as np

def test_height_k_value_optimization():
    """Test k-value optimization for height analysis."""
    print("🎯 Height Analysis K-Value Optimization Test")
    print("=" * 50)
    
    try:
        from .k_value_optimizer import KValueOptimizer
        
        # Create optimizer specifically for height analysis
        print("Creating k-value optimizer for height analysis...")
        optimizer = KValueOptimizer(
            analysis_type='height',
            n_trials=12,  # Reasonable number of trials
            k_range=(8, 30)  # Height analysis typically uses smaller k-values
        )
        
        # Run optimization
        print("Running k-value optimization for height analysis...")
        start_time = time.time()
        
        results = optimizer.optimize(
            gro_file="cases/lnb.gro",
            xtc_file="cases/md.xtc",
            residues={'DPPC': ['PO4'], 'CHOL': ['ROH']}
        )
        
        optimization_time = time.time() - start_time
        
        # Display results
        print(f"\n✅ Height analysis k-value optimization completed!")
        print(f"   Optimization time: {optimization_time:.2f} seconds")
        print(f"   Optimal k-value: {results['best_k']}")
        print(f"   Best score: {results['best_score']:.4f}")
        print(f"   K-values tested: {len(results['k_values_tested'])}")
        print(f"   Tested k-values: {sorted(results['k_values_tested'])}")
        
        # Get system features
        features = results['system_features']
        print(f"\n📊 System Features:")
        for key, value in features.items():
            print(f"   {key}: {value}")
        
        # Get recommendations
        recommendations = optimizer.get_recommendations(results)
        print(f"\n💡 Recommendations:")
        print(f"   Confidence: {recommendations['confidence']}")
        print(f"   Reasoning: {recommendations['reasoning']}")
        print(f"   Alternative k-values: {recommendations['alternative_k_values']}")
        
        # Test the optimal k-value manually
        print(f"\n🧪 Testing optimal k-value manually...")
        test_optimal_k(results['best_k'])
        
        # Test a few other k-values for comparison
        print(f"\n📈 Comparing with other k-values...")
        compare_k_values([results['best_k'] - 2, results['best_k'], results['best_k'] + 2])
        
        return results
        
    except Exception as e:
        print(f"❌ Error during height k-value optimization: {e}")
        import traceback
        traceback.print_exc()
        return None

def test_optimal_k(k_value):
    """Test the optimal k-value manually."""
    try:
        from analysis.height import Height
        
        print(f"   Testing k = {k_value}...")
        start_time = time.time()
        
        # Create height analyzer with optimal k-value
        import MDAnalysis as mda
        universe = mda.Universe("cases/lnb.gro", "cases/md.xtc")
        analyzer = Height(
            universe=universe,
            residuesGroup={'DPPC': [['PO4']], 'CHOL': [['ROH']]},
            k=k_value
        )
        
        # Run analysis
        results = analyzer.run()
        computation_time = time.time() - start_time
        
        print(f"     Computation time: {computation_time:.3f}s")
        print(f"     Results obtained: {len(results) if hasattr(results, '__len__') else 'N/A'}")
        
        return computation_time, results
        
    except Exception as e:
        print(f"     Error testing k = {k_value}: {e}")
        return None, None

def compare_k_values(k_values):
    """Compare different k-values."""
    results_comparison = {}
    
    for k in k_values:
        if k > 0:  # Only test positive k-values
            print(f"   Comparing k = {k}...")
            computation_time, results = test_optimal_k(k)
            if computation_time is not None:
                results_comparison[k] = {
                    'computation_time': computation_time,
                    'results_count': len(results) if hasattr(results, '__len__') else 0
                }
    
    # Display comparison
    if results_comparison:
        print(f"\n📊 K-Value Comparison:")
        print(f"   {'K-Value':<8} {'Time (s)':<10} {'Results':<10}")
        print(f"   {'-' * 30}")
        for k, data in sorted(results_comparison.items()):
            print(f"   {k:<8} {data['computation_time']:<10.3f} {data['results_count']:<10}")

def test_height_analysis_with_different_k():
    """Test height analysis with different k-values manually."""
    print(f"\n🔬 Manual Height Analysis Test")
    print("=" * 40)
    
    k_values_to_test = [10, 15, 20, 25]
    results_summary = {}
    
    for k in k_values_to_test:
        print(f"\nTesting height analysis with k = {k}...")
        
        try:
            from analysis.height import Height
            
            start_time = time.time()
            
            # Create analyzer
            import MDAnalysis as mda
            universe = mda.Universe("cases/lnb.gro", "cases/md.xtc")
            analyzer = Height(
                universe=universe,
                residuesGroup={'DPPC': [['PO4']], 'CHOL': [['ROH']]},
                k=k
            )
            
            # Run analysis
            results = analyzer.run()
            computation_time = time.time() - start_time
            
            results_summary[k] = {
                'computation_time': computation_time,
                'success': True,
                'results_count': len(results) if hasattr(results, '__len__') else 0
            }
            
            print(f"   ✅ Success! Time: {computation_time:.3f}s, Results: {results_summary[k]['results_count']}")
            
        except Exception as e:
            results_summary[k] = {
                'computation_time': None,
                'success': False,
                'error': str(e)
            }
            print(f"   ❌ Failed: {e}")
    
    # Display summary
    print(f"\n📋 Manual Test Summary:")
    print(f"   {'K-Value':<8} {'Status':<8} {'Time (s)':<10} {'Results':<10}")
    print(f"   {'-' * 40}")
    for k, data in sorted(results_summary.items()):
        status = "✅ Success" if data['success'] else "❌ Failed"
        time_str = f"{data['computation_time']:.3f}" if data['computation_time'] else "N/A"
        results_str = str(data['results_count']) if data['success'] else "N/A"
        print(f"   {k:<8} {status:<8} {time_str:<10} {results_str:<10}")

def main():
    """Main test function."""
    print("🚀 Height Analysis K-Value Optimization Test")
    print("=" * 60)
    
    # Test 1: K-value optimization
    print("Test 1: Automated K-Value Optimization")
    print("-" * 40)
    optimization_results = test_height_k_value_optimization()
    
    # Test 2: Manual comparison
    print("\nTest 2: Manual K-Value Comparison")
    print("-" * 40)
    test_height_analysis_with_different_k()
    
    # Summary
    print(f"\n🎉 Test Summary:")
    print(f"=" * 40)
    if optimization_results:
        print(f"✅ K-value optimization completed successfully!")
        print(f"   Recommended k-value: {optimization_results['best_k']}")
        print(f"   Optimization score: {optimization_results['best_score']:.4f}")
        print(f"   Total trials: {len(optimization_results['k_values_tested'])}")
    else:
        print(f"❌ K-value optimization failed!")
    
    print(f"\n💡 Next steps:")
    print(f"   1. Use the recommended k-value for your height analysis")
    print(f"   2. Run: python analysis/height.py --k-value {optimization_results['best_k'] if optimization_results else 'XX'} --gro-file cases/lnb.gro --xtc-file cases/md.xtc")
    print(f"   3. Compare results with other k-values if needed")

if __name__ == "__main__":
    main()
