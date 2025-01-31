#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 [-d] [-v] <target_string> <replacement_string>"
    echo "  -d  Dry run: show what would be changed without making any changes"
    echo "  -v  Verbose: output detailed information about the changes"
    exit 1
}

# Parse command line options
dry_run=false
verbose=false
while getopts "dv" opt; do
    case $opt in
        d) dry_run=true ;;
        v) verbose=true ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))

# Check if the correct number of arguments is provided
if [ $# -ne 2 ]; then
    usage
fi

target_string=$1
replacement_string=$2

# Function to count occurrences of a pattern in a file
count_occurrences() {
    grep -o "$1" "$2" | wc -l
}

# Loop through all files in the directory and subdirectories, excluding hidden files and folders
find . -type f ! -path '*/\.*' | while read -r file; do
    # Count occurrences before making changes
    count_target=$(count_occurrences "$target_string" "$file")

    # Perform replacements using sed if not a dry run
    if [ "$dry_run" = false ]; then
        source_file=/tmp/"$file_"tmp
        cp -f $file $source_file;
        outfile=$file
        sed -i "s/$target_string/$replacement_string/g" "$file"

    else
        source_file=$file;
        outfile=/tmp/"$file"_tmp
        sed "s/$target_string/$replacement_string/g" "$file" > $outfile
    fi

    # Count occurrences after making changes (or simulate for dry run)
    new_count_target=$(count_occurrences "$replacement_string" "$outfile")
 

    # Output the file name and number of changes if any changes were made (or simulated)
    if [ $count_target -gt 0 ]; then
        echo "Modified file: $file"
        echo "Changes made:"
        echo "  $target_string -> $replacement_string: $changes_target"

        # Output detailed information if verbose mode is enabled
        if [ "$verbose" = true & ]; then
            echo "Before:"
            grep --color=always -n "$target_string" "$source_file"
            echo "After:"
            grep --color=always -n "$replacement_string" "$outfile"
        fi
    fi
done

# If dry run, inform the user that no changes were actually made
if [ "$dry_run" = true ]; then
    echo "Dry run completed. No changes were actually made."
fi
