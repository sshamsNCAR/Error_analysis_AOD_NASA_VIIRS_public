import os
import pandas as pd
from datetime import datetime

# Define paths
input_dir = "/AOD_Level15_All_Points_V3/AOD/AOD15/ALL_POINTS/"
#level 2 data
input_dir = "/AOD_Level2_All_Points_V3/AOD20/ALL_POINTS/"
output_dir = "/AOD_Level2_All_Points_V3/formatted/"
# level = "lev15"  # Change to "1.5" for Level 1.5 data
level = "lev20"  # Change to "2" for Level 2 data

AERONET_WAVELENGTHS = [1.640, 1.020, 0.870, 0.865, 0.779, 0.675, 0.667, 0.620,
                       0.560, 0.555, 0.551, 0.532, 0.531, 0.510, 0.500, 0.490,
                       0.443, 0.440, 0.412, 0.400, 0.380, 0.340]

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

def convert_to_julian_day(date):
    """Convert a date (YYYY-MM-DD) to Julian day."""
    dt = datetime.strptime(date, "%Y-%m-%d")
    year = dt.year
    julian_day = dt.timetuple().tm_yday
    return year, julian_day

def process_aeronet_file(file_path):
    """Process a single .lev15 file and generate formatted .dat files."""
    # Extract station name from the file name
    file_name = os.path.basename(file_path)
    station_name = "_".join(file_name.split("_")[2:]).replace("."+level, "")
    

    # Read the .lev15 or lev 2 file
    try:
        data = pd.read_csv(file_path, skiprows=6, delimiter=",", engine="python",  encoding="ISO-8859-1")
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return

     # Ensure the file has the required columns
    required_columns = ["Date(dd:mm:yyyy)", "Time(hh:mm:ss)","Day_of_Year(Fraction)"] + [f"AOD_{int(w*1000)}nm" for w in AERONET_WAVELENGTHS]+["440-870_Angstrom_Exponent"]
    if not all(col in data.columns for col in required_columns):
        print(f"File {file_path} is missing required columns.")
        return

    # Process each row and write to the appropriate .dat file
    for _, row in data.iterrows():
        try:
            # Parse date and timeq
            date_str = row["Date(dd:mm:yyyy)"]
            time_str = row["Time(hh:mm:ss)"]
            day_of_year_fraction = row["Day_of_Year(Fraction)"]
            Angstrom_exponent = row["440-870_Angstrom_Exponent"]
            
             # Skip rows with missing AOD values for any wavelength
            if any(pd.isna(row[f"AOD_{int(w*1000)}nm"]) for w in AERONET_WAVELENGTHS):
                continue

            # Convert date to Julian day
            date_obj = datetime.strptime(date_str, "%d:%m:%Y")
            year, julian_day = date_obj.year, date_obj.timetuple().tm_yday

            # Create output directory for the year and Julian day
            year_dir = os.path.join(output_dir, str(year))
            julian_day_dir = os.path.join(year_dir, f"{julian_day:03d}")
            os.makedirs(julian_day_dir, exist_ok=True)

            # Write to the station's .dat file
            output_file = os.path.join(julian_day_dir, f"{station_name}.dat")
            with open(output_file, "a") as f:
                # Write all wavelengths and their AOD values
                aod_values = [row[f"AOD_{int(w*1000)}nm"] for w in AERONET_WAVELENGTHS]
                aod_values_str = ",".join(map(str, aod_values))
                f.write(f"{date_str},{time_str},{day_of_year_fraction:.6f},{aod_values_str},{Angstrom_exponent}\n")
        except Exception as e:
            print(f"Error processing row in file {file_path}: {e}")

def main():
    # Process all aeronet data (.lev15 or .lev20) files in the input directory
    all_files = []
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith("." + level):
                all_files.append(os.path.join(root, file))
                #process_aeronet_file(file_path)
     # Sort files for consistent processing
    all_files.sort()

    # Process files in batches
    #batch_size = 50
    start_index = 800  # Start from file 50q
    end_index = 1800   # End at file 100

    # Process the first batch (files 50 to 100)
    for file_path in all_files[start_index:end_index]:
        process_aeronet_file(file_path)
        print(file_path)

    # # Process the remaining files in batches
    # for i in range(end_index, len(all_files), batch_size):
    #     batch_files = all_files[i:i + batch_size]
    #     for file_path in batch_files

    #         process_aeronet_file(file_path)
    #         print(file_path)

if __name__ == "__main__":
    main()