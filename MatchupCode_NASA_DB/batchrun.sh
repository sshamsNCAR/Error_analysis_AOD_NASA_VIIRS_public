#!/bin/bash
# example usage: ./batchrun.sh 2021 21 2021 22 > ../Bin/DT_batchrun.log 2>&1
# filepath: /glade/u/home/sshams/code/AFRICA_SERVIR/Error_analysis_AOD_NASA/MatchupCode_NASA_DT/batchrun.sh

# Paths
AFRICA_LIST="/glade/campaign/acom/acom-da/SERVIR/ind-obs/aeronet/aeronet_stations_africa_2021_2025.txt"
#DATA_DIR="/glade/campaign/acom/acom-da/SERVIR/ind-obs/aeronet/AOD_Level15_All_Points_V3/for/"
DATA_DIR="/glade/campaign/acom/acom-da/SERVIR/ind-obs/aeronet/AOD_Level2_All_Points_V3/formatted/"
#OUTPUT_DIR="/glade/campaign/acom/acom-da/SERVIR/ind-obs/aeronet/AOD_Level15_All_Points_V3/for/"
OUTPUT_DIR="/glade/campaign/acom/acom-da/SERVIR/ind-obs/aeronet/AOD_Level2_All_Points_V3/formatted/"

MATCH_code_dir="/glade/u/home/sshams/code/AFRICA_SERVIR/Error_analysis_AOD_NASA/Bin/match_viirs_db_n20"
# exec > /glade/u/home/sshams/code/AFRICA_SERVIR/Error_analysis_AOD_NASA/Bin/match_viirs_dt_n20/DB_batchrun.log 2>&1



# Read input arguments for start and end year/day
START_YEAR=$1
START_DAY=$2
END_YEAR=$3
END_DAY=$4

if [[ -z "$START_YEAR" || -z "$START_DAY" || -z "$END_YEAR" || -z "$END_DAY" ]]; then
    echo "Usage: $0 <start_year> <start_day> <end_year> <end_day>"
    exit 1
fi

# Read the Africa station list into an array
declare -A africa_stations
while IFS=',' read -r siteid longitude latitude elevation; do
    africa_stations["$siteid"]=1
done < <(tail -n +2 "$AFRICA_LIST")  # Skip the header line
echo "Africa station list loaded with ${#africa_stations[@]} stations."

# Iterate through the specified range of years and days
current_year=$START_YEAR
current_day=$START_DAY

while [[ $current_year -lt $END_YEAR || ($current_year -eq $END_YEAR && $current_day -le $END_DAY) ]]; do
    # Format the day as a 3-digit number (e.g., 001, 002, ..., 366)
    day_formatted=$(printf "%03d" $current_day)

    # Path to the current day's directory
    day_dir="${DATA_DIR}${current_year}/${day_formatted}/"

    if [[ -d "$day_dir" ]]; then
        echo "Processing year $current_year, day $day_formatted..."

        # Output file for the current day
        output_file="${OUTPUT_DIR}aeronet_locations.txt"

        # Initialize the output file with a header
        echo "siteid,longitude,latitude,elevation" > "$output_file"

        # Check for station files in the current day folder
        if compgen -G "$day_dir"*.dat > /dev/null; then
            for station_file in "$day_dir"*.dat; do
                station_name=$(basename "$station_file" .dat)

                # If the station is in the Africa list, add it to the output file
                if [[ -n "${africa_stations[$station_name]}" ]]; then
                    # Extract station details from the Africa list
                    grep "^$station_name," "$AFRICA_LIST" >> "$output_file"
                fi
            done
        else
            echo "    No station files found in $day_dir. Skipping..."
        fi

        # Run the match_viirs_dt_n20 command for the current day
        # month=$(printf "%02d" $(( (current_day - 1) / 30 + 1 )))  # Approximate month
        # day_of_month=$(printf "%02d" $(( (current_day - 1) % 30 + 1 )))  # Approximate day
        # Convert DOY (day of year) to a valid date using the `date` command
        date_str=$(date -d "$current_year-01-01 +$((current_day - 1)) days" +"%Y %m %d")
        read year month day_of_month <<< "$date_str"
        echo "Running match_viirs_db_n20 for $current_year-$month-$day_of_month..."
        $MATCH_code_dir $current_year $month $day_of_month
        # Wait for the command to finish before proceeding
        wait
    else
        echo "Directory not found for year $current_year, day $day_formatted. Skipping..."
    fi

    # Increment the day
    current_day=$((current_day + 1))

    # Check if we need to move to the next year
    if [[ $current_day -gt 366 ]]; then
        current_day=1
        current_year=$((current_year + 1))
    fi
done
