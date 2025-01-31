#!/bin/bash

# Script Name: rename_files.sh
# Description: This script renames files in a directory and its subdirectories whose names contain a target string, replacing it with a replacement string. It supports dry run, verbosity, using 'git mv' if the file is under version control, and using extended regex.
# Creation Date: 2025-01-25
# Version: 1.1

# Function to display usage information
usage() {
    echo "Usage: $0 [-d] [-v] [-g] [-e] <target_string> <replacement_string>"
    echo "  -d  Dry run: show what would be changed without making any changes"
    echo "  -v  Verbose: output detailed information about the changes"
    echo "  -g  Use 'git mv' instead of 'mv' to rename files if they are under version control"
    echo "  -e  Use extended regex for matching the target string"
    exit 1
}

# Parse command line options
dry_run=false
verbose=false
use_git_mv=false
extended_regex=false
while getopts "dvge" opt; do
    case $opt in
        d) dry_run=true ;;
        v) verbose=true ;;
        g) use_git_mv=true ;;
        e) extended_regex=true ;;
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

# Function to check if a file is under version control (git)
is_git_controlled() {
    git ls-files --error-unmatch "$1" > /dev/null 2>&1
}

# Loop through all files in the directory and subdirectories, excluding hidden files and folders
find . -type f ! -path '*/\.*' | while read -r file; do
    # Check if the file name contains the target string (using extended regex if specified)
    if [ "$extended_regex" = true ]; then
        if [[ "$file" =~ $target_string ]]; then
            new_file=$(echo "$file" | sed -E "s/$target_string/$replacement_string/g")
        else
            continue
        fi
    else
        if [[ "$file" == *"$target_string"* ]]; then
            new_file=$(echo "$file" | sed "s/$target_string/$replacement_string/g")
        else
            continue
        fi
    fi

    # Output the changes if verbose mode is enabled or if it's a dry run
    if [ "$verbose" = true ] || [ "$dry_run" = true ]; then
        echo "Renaming: $file -> $new_file"
    fi

    # Perform the renaming if not a dry run
    if [ "$dry_run" = false ]; then
        if [ "$use_git_mv" = true ] && is_git_controlled "$file"; then
            git mv "$file" "$new_file"
        else
            mv "$file" "$new_file"
        fi
    fi
done

# If dry run, inform the user that no changes were actually made
if [ "$dry_run" = true ]; then
    echo "Dry run completed. No changes were actually made."
fi
