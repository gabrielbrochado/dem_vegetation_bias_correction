
!pip install pingouin
import pandas as pd
import pingouin as pg

# Read the CSV file into a pandas DataFrame
base_path = 'C:/Users/gabro/Downloads' # <- where suplementary materials was unziped
file_path = base_path + '/Suplementary_Materials_Brochado_Renno_2024/stacked_dataframe_r_class_p.csv'
df = pd.read_csv(file_path)

# Define the columns for the Wilcoxon test
columns_to_test = [
    ('glo_corrected_area', 'glo_area'),
    ('glo_corrected_area', 'fabdem_area'),
    ('fabdem_area', 'glo_area')
]

# Define the area and window combinations
areas = [1, 2, 3, 4]
radii = [1000, 2000, 3000]

# Number of total samples to select
total_samples = 50

# Initialize a list to store results and to accumulate the selected DataFrame
results = []
compiled_selected_df = pd.DataFrame()  # Empty DataFrame to store all selected samples

# Iterate over each area and radius combination
for area in areas:
    for radius in radii:
        # Filter the DataFrame for the current area and radius and make a copy
        filtered_df = df[(df['area'] == area) & (df['radius'] == radius)].copy()

        # Skip if the filtered DataFrame is empty
        if filtered_df.empty:
            results.append({
                'area': area,
                'radius': radius,
                'corrected x glo': None,
                'corrected x fabdem': None,
                'fabdem x glo': None,
                'total_flowpaths': 0,
                'class_0_count': 0,
                'class_1_count': 0,
                'samples_taken': 0,
                'samples_class_0': 0,
                'samples_class_1': 0,
                'e_CxG_samples': 0,
                'e_CxF_samples': 0,
                'e_FxG_samples': 0
            })
            continue

        # Create 'min_area' field with the minimum among the specified fields
        filtered_df['min_area'] = filtered_df[['glo_corrected_area', 'glo_area', 'fabdem_area']].min(axis=1)

        # Determine the total number of flowpaths and count for each class
        class_0_df = filtered_df[filtered_df['class'] == 0].copy()
        class_1_df = filtered_df[filtered_df['class'] == 1].copy()

        class_0_count = len(class_0_df)
        class_1_count = len(class_1_df)

        # Sort the DataFrames by 'min_area' in ascending order
        class_0_df = class_0_df.sort_values(by='min_area').copy()
        class_1_df = class_1_df.sort_values(by='min_area').copy()

        # Calculate the number of samples to take from each class proportionally
        if class_0_count + class_1_count <= total_samples:
            # If total available is less than or equal to the number we want, take all
            selected_class_0_df = class_0_df
            selected_class_1_df = class_1_df
        else:
            # Calculate proportional sample sizes based on the counts
            if class_0_count > 0 and class_1_count > 0:
                samples_class_0 = round((class_0_count / (class_0_count + class_1_count)) * total_samples)
                samples_class_1 = total_samples - samples_class_0
            elif class_0_count > 0:
                samples_class_0 = total_samples
                samples_class_1 = 0
            else:
                samples_class_0 = 0
                samples_class_1 = total_samples

            # Select samples based on calculated sizes
            selected_class_0_df = class_0_df.iloc[:samples_class_0].copy()
            selected_class_1_df = class_1_df.iloc[:samples_class_1].copy()

        # Combine the sampled data
        selected_df = pd.concat([selected_class_0_df, selected_class_1_df])

        # Add area and radius columns to selected_df for identification
        selected_df['area'] = area
        selected_df['radius'] = radius

        # Append the current selected_df to the compiled_selected_df
        compiled_selected_df = pd.concat([compiled_selected_df, selected_df])

        # Initialize a dictionary to store p-values for this area and radius
        p_values = {}

        # Initialize counts for effective samples (excluding ties)
        e_CxG_samples = 0
        e_CxF_samples = 0
        e_FxG_samples = 0

        # Perform the Wilcoxon test for each pair of columns on the selected rows
        for col1, col2 in columns_to_test:
            # Convert columns to numeric and drop rows with non-numeric values
            selected_df[col1] = pd.to_numeric(selected_df[col1], errors='coerce')
            selected_df[col2] = pd.to_numeric(selected_df[col2], errors='coerce')
            selected_df = selected_df.dropna(subset=[col1, col2])

            # Count effective samples excluding ties
            if col1 == 'glo_corrected_area' and col2 == 'glo_area':
                e_CxG_samples = len(selected_df[selected_df[col1] != selected_df[col2]])
            elif col1 == 'glo_corrected_area' and col2 == 'fabdem_area':
                e_CxF_samples = len(selected_df[selected_df[col1] != selected_df[col2]])
            elif col1 == 'fabdem_area' and col2 == 'glo_area':
                e_FxG_samples = len(selected_df[selected_df[col1] != selected_df[col2]])

            # Skip if there are not enough samples to perform the test
            if len(selected_df) < 2:
                p_values[f'{col1} x {col2}'] = None
                continue

            # Perform Wilcoxon test using Pingouin and get p-value
            wilcoxon_result = pg.wilcoxon(selected_df[col1], selected_df[col2], alternative='less')
            p_value = wilcoxon_result['p-val'].values[0]

            # Check if the p-value is greater than 0.5
            if p_value > 0.5:
                p_value = 1 - p_value
                p_value = f"{round(p_value, 6)}*"
            else:
                p_value = round(p_value, 6)

            p_values[f'{col1} x {col2}'] = p_value

        # Store the results for the current area and radius
        results.append({
            'area': area,
            'radius': radius,
            'corrected x glo': p_values['glo_corrected_area x glo_area'],
            'corrected x fabdem': p_values['glo_corrected_area x fabdem_area'],
            'fabdem x glo': p_values['fabdem_area x glo_area'],
            'total_flowpaths': len(filtered_df),
            'class_0_count': class_0_count,
            'class_1_count': class_1_count,
            'samples_taken': len(selected_df),
            'samples_class_0': len(selected_class_0_df),
            'samples_class_1': len(selected_class_1_df),
            'e_CxG_samples': e_CxG_samples,
            'e_CxF_samples': e_CxF_samples,
            'e_FxG_samples': e_FxG_samples
        })

# Convert the results list to a DataFrame for display
results_df = pd.DataFrame(results)

# Display the compiled selected DataFrame and the results DataFrame
print("Compiled Selected DataFrame:")

compiled_selected_df.to_csv(base_path + '/Suplementary_Materials_Brochado_Renno_2024/selected_stacked_dataframe_r_class_p.csv', index=False)
print("\nP-values DataFrame:")
results_df
