#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Command-line figure generation tool - Automatically detect available plot types from CSV files and generate figures

Usage:
    python figure/make_figure.py -f area.csv -o ./cases/figure.png
"""

import os
import sys
import argparse
import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端，适合命令行

import matplotlib.pyplot as plt
import pandas as pd

# 添加项目根目录到路径
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from figure.figure import read_excel, LipidsFigure, BubbleFigure, TYPE
from figure.density_figure import DensityFigure


def get_available_plot_types(description, data_type):
    """
    Determine available plot types based on description and data_type
    
    Args:
        description: Description from CSV file
        data_type: Type value from TYPE dictionary
    
    Returns:
        list: List of available plot types ['line', 'bar', 'scatter', 'heatmap']
    """
    available_types = []
    
    # Check if it's a Density type based on description keywords
    # This handles cases where description is not in TYPE dictionary
    if any(keyword in description for keyword in ['DensityTime', 'Density Radius', 'Density With Times', 
                                                  'Density With Time', 'Density With Radius', 'Multi-Radius Density',
                                                  'Multi Radius Density', 'Density With']):
        # DensityFigure: supports Line, Heatmap
        available_types = ['line', 'heatmap']
        return available_types
    
    # Determine based on data_type
    if data_type == 0:
        # LipidsFigure: supports Line, Bar, Scatter
        available_types = ['line', 'bar', 'scatter']
    elif data_type == 1:
        # BubbleFigure: supports Line, Bar
        available_types = ['line', 'bar']
    elif data_type == 3:
        # DensityFigure: supports Line, Heatmap
        available_types = ['line', 'heatmap']
    elif data_type == 2:
        # Other types, default to Line, Bar
        available_types = ['line', 'bar']
    else:
        # Unknown type, default to Line, Bar
        available_types = ['line', 'bar']
    
    return available_types


def create_figure_instance(description, excel_data, figure_settings):
    """
    Create corresponding figure instance based on description
    
    Args:
        description: Description from CSV file
        excel_data: pandas DataFrame
        figure_settings: Figure settings dictionary
    
    Returns:
        Figure instance
    """
    # Check if it's a Density type
    if any(keyword in description for keyword in ['DensityTime', 'Density Radius', 'Density With Times', 
                                                  'Density With Time', 'Density With Radius', 'Multi-Radius Density']):
        return DensityFigure(description, excel_data, figure_settings)
    
    # Determine based on TYPE dictionary
    if description in TYPE:
        data_type = TYPE[description]
        
        if data_type == 0:
            return LipidsFigure(description, excel_data, figure_settings)
        elif data_type == 1:
            return BubbleFigure(description, excel_data, figure_settings)
        elif data_type == 3:
            return DensityFigure(description, excel_data, figure_settings)
        else:
            # Default to LipidsFigure
            return LipidsFigure(description, excel_data, figure_settings)
    else:
        # Infer from data structure
        # If contains Resname column, use LipidsFigure
        if 'Resname' in excel_data.columns:
            return LipidsFigure(description, excel_data, figure_settings)
        # Otherwise use BubbleFigure
        else:
            return BubbleFigure(description, excel_data, figure_settings)


def interactive_plot_selection(available_types):
    """
    Interactive plot type selection
    
    Args:
        available_types: List of available plot types
    
    Returns:
        str: Selected plot type or 'all'
    """
    print("\nAvailable plot types:")
    for i, plot_type in enumerate(available_types, 1):
        type_names = {
            'line': 'Line Chart',
            'bar': 'Bar Chart',
            'scatter': 'Scatter Chart',
            'heatmap': 'Heatmap Chart'
        }
        print(f"  {i}. {type_names.get(plot_type, plot_type)}")
    print(f"  {len(available_types) + 1}. Generate All")
    
    while True:
        try:
            choice = input(f"\nPlease select plot type (1-{len(available_types) + 1}): ").strip()
            choice_num = int(choice)
            
            if 1 <= choice_num <= len(available_types):
                return available_types[choice_num - 1]
            elif choice_num == len(available_types) + 1:
                return 'all'
            else:
                print(f"Invalid choice, please enter a number between 1 and {len(available_types) + 1}")
        except ValueError:
            print("Please enter a valid number")
        except KeyboardInterrupt:
            print("\n\nOperation cancelled")
            sys.exit(0)


def generate_plot(figure_instance, plot_type, output_path, base_name):
    """
    Generate and save figures
    
    Args:
        figure_instance: Figure instance
        plot_type: Plot type ('line', 'bar', 'scatter', 'heatmap', 'all')
        output_path: Output directory or file path
        base_name: Base filename (without extension)
    
    Returns:
        list: List of generated file paths
    """
    generated_files = []
    
    # Ensure output directory exists
    output_dir = os.path.dirname(output_path) if os.path.dirname(output_path) else '.'
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    # Determine plot types to generate
    plot_types_to_generate = []
    if plot_type == 'all':
        # Get all methods supported by the instance
        if hasattr(figure_instance, 'plot_line'):
            plot_types_to_generate.append('line')
        if hasattr(figure_instance, 'plot_bar'):
            plot_types_to_generate.append('bar')
        if hasattr(figure_instance, 'plot_scatter'):
            plot_types_to_generate.append('scatter')
        if hasattr(figure_instance, 'plot_heatmap'):
            plot_types_to_generate.append('heatmap')
    else:
        plot_types_to_generate = [plot_type]
    
    # Save original plt.show function (if exists)
    original_show = plt.show
    
    # Temporarily replace plt.show to avoid popup windows in command-line mode
    def noop_show(*args, **kwargs):
        pass
    
    # Generate plots
    for pt in plot_types_to_generate:
        plt.figure(figsize=(10, 6))
        
        try:
            # Temporarily replace plt.show
            plt.show = noop_show
            
            # Call corresponding plot method
            if pt == 'line':
                figure_instance.plot_line()
            elif pt == 'bar':
                figure_instance.plot_bar()
            elif pt == 'scatter':
                if hasattr(figure_instance, 'plot_scatter'):
                    figure_instance.plot_scatter()
                else:
                    print(f"Warning: Plot type '{pt}' is not supported, skipping")
                    plt.close()
                    continue
            elif pt == 'heatmap':
                if hasattr(figure_instance, 'plot_heatmap'):
                    figure_instance.plot_heatmap()
                else:
                    print(f"Warning: Plot type '{pt}' is not supported, skipping")
                    plt.close()
                    continue
            else:
                print(f"Warning: Unknown plot type '{pt}', skipping")
                plt.close()
                continue
            
            # Restore plt.show
            plt.show = original_show
            
            # Generate output filename
            if plot_type == 'all':
                # If 'all', generate files with suffix for each type
                output_file = os.path.join(output_dir, f"{base_name}_{pt}.png")
            else:
                # If only one type, use original output path
                if len(plot_types_to_generate) == 1:
                    # Ensure output path has .png extension
                    if not output_path.endswith('.png'):
                        output_file = output_path + '.png'
                    else:
                        output_file = output_path
                else:
                    # If multiple types generated (should not happen), add suffix
                    output_file = os.path.join(output_dir, f"{base_name}_{pt}.png")
            
            # Save image
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            generated_files.append(output_file)
            print(f"✓ Figure saved: {output_file}")
            
        except Exception as e:
            # Restore plt.show
            plt.show = original_show
            print(f"✗ Error generating {pt} plot: {e}")
            import traceback
            traceback.print_exc()
            plt.close()
            continue
    
    # Ensure plt.show is restored
    plt.show = original_show
    
    return generated_files


def parse_args():
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser(
        description="Command-line figure generation tool - Automatically detect available plot types from CSV files and generate figures",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python figure/make_figure.py -f area.csv -o ./cases/figure.png
  python figure/make_figure.py -f anisotropy_results.csv -o ./cases/result.png --plot-type line
  python figure/make_figure.py -f height_results.csv -o ./cases/plots --plot-type all
        """
    )
    
    parser.add_argument(
        '-f', '--file',
        type=str,
        required=True,
        help='Input CSV file path'
    )
    
    parser.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help='Output image file path (will automatically add suffix based on plot type)'
    )
    
    parser.add_argument(
        '--plot-type',
        type=str,
        choices=['line', 'bar', 'scatter', 'heatmap', 'all'],
        default=None,
        help='Plot type. If not specified, an interactive menu will be displayed for user selection'
    )
    
    parser.add_argument(
        '--no-interactive',
        action='store_true',
        help='Non-interactive mode. If --plot-type is not specified, will use default type (line)'
    )
    
    return parser.parse_args()


def main():
    """Main function"""
    args = parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.file):
        print(f"Error: Input file does not exist: {args.file}")
        sys.exit(1)
    
    # Read CSV file
    print(f"Reading file: {args.file}")
    try:
        description, excel_data, time_unit = read_excel(args.file)
        print(f"Data description: {description}")
        print(f"Data shape: {excel_data.shape}")
        print(f"Time unit: {time_unit}")
    except Exception as e:
        print(f"Error: Failed to read CSV file: {e}")
        sys.exit(1)
    
    # Determine data type and available plot types
    if description in TYPE:
        data_type = TYPE[description]
    else:
        # Check if it's a Density type based on description keywords
        if any(keyword in description for keyword in ['DensityTime', 'Density Radius', 'Density With Times', 
                                                      'Density With Time', 'Density With Radius', 'Multi-Radius Density',
                                                      'Multi Radius Density', 'Density With']):
            data_type = 3  # DensityFigure
        # Infer from data structure
        elif 'Resname' in excel_data.columns:
            data_type = 0  # LipidsFigure
        else:
            data_type = 1  # BubbleFigure
    
    available_types = get_available_plot_types(description, data_type)
    
    if not available_types:
        print("Error: Unable to determine available plot types")
        sys.exit(1)
    
    # Determine plot type
    if args.plot_type:
        plot_type = args.plot_type
        if plot_type not in available_types and plot_type != 'all':
            print(f"Error: Specified plot type '{plot_type}' is not available")
            print(f"Available types: {', '.join(available_types)}")
            sys.exit(1)
    elif args.no_interactive:
        plot_type = available_types[0]  # Default to first available type
        print(f"Non-interactive mode, using default type: {plot_type}")
    else:
        plot_type = interactive_plot_selection(available_types)
    
    # Create figure instance
    figure_settings = {
        'x_title': 'Time (ns)',
        'y_title': description,
        'axis_text': 12
    }
    
    try:
        figure_instance = create_figure_instance(description, excel_data, figure_settings)
    except Exception as e:
        print(f"Error: Failed to create figure instance: {e}")
        sys.exit(1)
    
    # Generate output filename (remove extension as base name)
    output_path_full = os.path.abspath(args.output)
    output_dir = os.path.dirname(output_path_full) if os.path.dirname(output_path_full) else '.'
    base_name = os.path.splitext(os.path.basename(output_path_full))[0]
    
    # If base_name is empty (user only specified directory), use input filename
    if not base_name:
        input_base = os.path.splitext(os.path.basename(args.file))[0]
        base_name = input_base
        output_path = os.path.join(output_dir, base_name) if output_dir else base_name
    else:
        output_path = os.path.join(output_dir, base_name) if output_dir else base_name
    
    # Generate plots
    print(f"\nGenerating {plot_type} plot(s)...")
    generated_files = generate_plot(figure_instance, plot_type, output_path, base_name)
    
    if generated_files:
        print(f"\n✓ Successfully generated {len(generated_files)} figure file(s):")
        for file_path in generated_files:
            print(f"  - {file_path}")
    else:
        print("\n✗ Failed to generate any figures")
        sys.exit(1)


if __name__ == '__main__':
    main()

