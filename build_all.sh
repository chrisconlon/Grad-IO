#!/bin/bash

# Build script for all LaTeX files in Grad-IO project
# Uses XeLaTeX with latexmk for proper handling of citations and cross-references

# Configuration
PROJECT_ROOT="$(pwd)"
BUILD_ENGINE="xelatex"
LOG_FILE="$PROJECT_ROOT/build_log.txt"
TEXLIVE_BIN="/Library/TeX/texbin/"

# Add TeX Live to PATH for tools like kpsewhich
export PATH="$TEXLIVE_BIN:$PATH"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1"
}

print_success() {
    echo -e "${GREEN}✓${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}⚠${NC} $1"
}

print_error() {
    echo -e "${RED}✗${NC} $1"
}

# Initialize log file
echo "Build started at $(date)" > "$LOG_FILE"

# Change to project directory
cd "$PROJECT_ROOT" || exit 1

print_status "Starting build process for all LaTeX files..."
print_status "Project root: $PROJECT_ROOT"
print_status "Build engine: $BUILD_ENGINE"
print_status "Log file: $LOG_FILE"

# Initialize counters
total_files=0
successful_builds=0
failed_builds=0

# Find all .tex files that are likely main documents (not includes)
# Exclude files that are typically included (like preamble, resources, etc.)
print_status "Scanning for LaTeX documents that use teaching_slides.sty..."

# Create array of .tex files to process
declare -a tex_files

# Find .tex files but exclude common include patterns
find "$PROJECT_ROOT" -name "*.tex" -type f | while read -r file; do
    # Skip files in resources directory
    if [[ "$file" == *"/resources/"* ]]; then
        continue
    fi
    
    # Skip files that look like includes (common patterns)
    basename_file=$(basename "$file")
    if [[ "$basename_file" == "preamble"* ]] || \
       [[ "$basename_file" == "header"* ]] || \
       [[ "$basename_file" == "footer"* ]] || \
       [[ "$basename_file" == "include"* ]]; then
        continue
    fi
    
    # Check if file contains \documentclass (indicating it's a main document)
    # AND uses teaching_slides package
    if grep -q "documentclass" "$file" && grep -q "teaching_slides" "$file"; then
        echo "$file"
    fi
done > /tmp/tex_files_list.txt

# Read the files into array using a more compatible method
while IFS= read -r line; do
    tex_files+=("$line")
done < /tmp/tex_files_list.txt

total_files=${#tex_files[@]}
print_status "Found $total_files LaTeX documents using teaching_slides.sty"

# Debug: show first few files found
if [ $total_files -gt 0 ]; then
    print_status "Sample files found:"
    for i in "${!tex_files[@]}"; do
        if [ $i -lt 3 ]; then
            rel_path=${tex_files[$i]#$PROJECT_ROOT/}
            echo "  - $rel_path"
        fi
    done
    if [ $total_files -gt 3 ]; then
        echo "  ... and $((total_files - 3)) more"
    fi
else
    print_error "No .tex files with \\documentclass found!"
    print_status "Debugging: Let's check what files exist..."
    find "$PROJECT_ROOT" -name "*.tex" -type f | head -5
    exit 1
fi

# Build each file
for tex_file in "${tex_files[@]}"; do
    # Get relative path for nicer output
    rel_path=${tex_file#$PROJECT_ROOT/}
    
    print_status "Building: $rel_path"
    
    # Get directory and filename
    tex_dir=$(dirname "$tex_file")
    tex_name=$(basename "$tex_file" .tex)
    
    # Change to the directory containing the .tex file
    cd "$tex_dir" || {
        print_error "Failed to change to directory: $tex_dir"
        ((failed_builds++))
        continue
    }
    
    # Build with latexmk using XeLaTeX
    if "$TEXLIVE_BIN/latexmk" -pdf -pdflatex="$TEXLIVE_BIN/$BUILD_ENGINE -interaction=nonstopmode" -quiet "$tex_name.tex" >> "$LOG_FILE" 2>&1; then
        print_success "Successfully built: $rel_path"
        ((successful_builds++))
    else
        print_error "Failed to build: $rel_path"
        echo "Error building $rel_path at $(date)" >> "$LOG_FILE"
        ((failed_builds++))
    fi
    
    # Return to project root
    cd "$PROJECT_ROOT" || exit 1
done

# Summary
echo ""
print_status "Build process completed!"
print_status "Total files processed: $total_files"
print_success "Successful builds: $successful_builds"
if [ $failed_builds -gt 0 ]; then
    print_error "Failed builds: $failed_builds"
else
    print_success "Failed builds: $failed_builds"
fi

# Show log file location
echo ""
print_status "Detailed log available at: $LOG_FILE"

# If there were failures, offer to show them
if [ $failed_builds -gt 0 ]; then
    echo ""
    print_warning "To see detailed error messages for failed builds:"
    echo "  tail -50 $LOG_FILE"
    echo ""
    print_warning "To build a specific file manually:"
    echo "  cd /path/to/tex/file"
    echo "  $TEXLIVE_BIN/latexmk -pdf -pdflatex='$TEXLIVE_BIN/$BUILD_ENGINE -interaction=nonstopmode' filename.tex"
fi

echo ""
echo "Build process finished at $(date)"