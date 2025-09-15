#!/usr/bin/env python3
"""
K-Value Optimization Example
===========================

This example shows how to automatically find the optimal k-value
for your molecular dynamics analysis.
"""

def main():
    """Main example function."""
    print("🎯 K-Value Optimization Example")
    print("=" * 40)
    
    try:
        from .k_value_optimizer import KValueOptimizer
        
        # Step 1: Create optimizer for your analysis type
        print("Step 1: Creating k-value optimizer...")
        optimizer = KValueOptimizer(
            analysis_type='area',  # Change this to your analysis type
            n_trials=15,           # Number of optimization trials
            k_range=(5, 40)        # Range of k-values to test
        )
        
        # Step 2: Run optimization
        print("Step 2: Running optimization...")
        results = optimizer.optimize(
            gro_file="cases/lnb.gro",      # Your GRO file
            xtc_file="cases/md.xtc",       # Your XTC file
            residues={'DPPC': ['PO4'], 'CHOL': ['ROH']}  # Your residues
        )
        
        # Step 3: Get results
        print("Step 3: Optimization completed!")
        print(f"   Optimal k-value: {results['best_k']}")
        print(f"   Best score: {results['best_score']:.4f}")
        print(f"   K-values tested: {len(results['k_values_tested'])}")
        
        # Step 4: Get recommendations
        recommendations = optimizer.get_recommendations(results)
        print(f"   Confidence: {recommendations['confidence']}")
        print(f"   Reasoning: {recommendations['reasoning']}")
        
        # Step 5: Use the optimal k-value
        print("Step 4: Use the optimal k-value in your analysis:")
        print(f"   python analysis/area.py --k-value {results['best_k']} --gro-file cases/lnb.gro --xtc-file cases/md.xtc")
        
        # Optional: Plot results
        print("Step 5: Plotting optimization results...")
        optimizer.plot_optimization_results(results, save_path="k_value_optimization.png")
        
        print("\n✅ K-value optimization completed successfully!")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        print("Make sure you have the required dependencies installed:")
        print("pip install scikit-learn matplotlib seaborn")

if __name__ == "__main__":
    main()
