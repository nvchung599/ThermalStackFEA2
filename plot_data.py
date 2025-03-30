import pandas as pd
import matplotlib.pyplot as plt
import sys

def plot_csv(file_path):
    print(f"Reading data from {file_path}...")
    
    # Read the CSV file
    try:
        # First check the structure with only a few rows
        df_preview = pd.read_csv(file_path, nrows=5)
        print("CSV structure preview:")
        print(df_preview.head())
        
        # Now read the full file
        df = pd.read_csv(file_path)
        print(f"Successfully loaded {len(df)} rows of data")
        
        # Basic info about the dataframe
        print("\nData columns:")
        for col in df.columns:
            print(f"- {col}")
        
        # Create a figure for dual Y-axis plot
        fig, ax1 = plt.subplots(figsize=(12, 8))
        
        # Define x-axis (time)
        x_column = 'time'
        
        # Set color cycle for temperature plots
        color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
        color_idx = 0
        
        # Plot temperature columns on left Y-axis (ax1)
        ax1.set_xlabel('Time')
        ax1.set_ylabel('Temperature (°C)')
        ax1.set_ylim(0, 100)  # Set temperature y-axis range to 0-100
        
        temp_columns = [col for col in df.columns if col.startswith('T_') and col != x_column]
        for col in temp_columns:
            color = color_cycle[color_idx % len(color_cycle)]
            color_idx += 1
            ax1.plot(df[x_column], df[col], label=col, linewidth=1.5, alpha=0.8)
        
        ax1.grid(True)
        
        # Create second Y-axis for power
        ax2 = ax1.twinx()
        
        # Plot power on right Y-axis (ax2)
        color = 'red'  # Distinct color for power
        power_col = 'P_resistor'
        ax2.set_ylabel('Power (W)', color=color)
        ax2.plot(df[x_column], df[power_col], color=color, linewidth=2, label=power_col)
        ax2.tick_params(axis='y', labelcolor=color)
        ax2.set_ylim(0, 1.1)  # Set power y-axis range to 0-1.1
        
        # Create combined legend
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left')
        
        plt.title(f'Temperature and Power vs Time from {file_path}')
        plt.tight_layout()
        
        # Save the plot
        output_file = "data_plot.png"
        plt.savefig(output_file)
        print(f"Plot saved to {output_file}")
        
        # Show the plot
        plt.show()
        
    except Exception as e:
        print(f"Error: {e}")
        return False
    
    return True

if __name__ == "__main__":
    file_path = "data.csv"
    
    # Use command line argument if provided
    if len(sys.argv) > 1:
        file_path = sys.argv[1]
    
    plot_csv(file_path)
