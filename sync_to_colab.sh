#!/bin/bash

# =============================================================================
# Sync Script: Local ↔ Colab
# =============================================================================
# This script helps sync your CFD Solver project between local VS Code and
# Google Colab remote instance
#
# Usage:
#   ./sync_to_colab.sh push    # Upload local changes to Colab
#   ./sync_to_colab.sh pull    # Download Colab changes to local
#   ./sync_to_colab.sh status  # Check sync status
#
# Prerequisites:
#   - SSH connection to Colab established
#   - SSH config entry for Colab host (see setup below)
# =============================================================================

# Configuration
LOCAL_PROJECT_DIR="/Users/rameshkolluru/My_Research/CFD_Solver_withCUDA"
REMOTE_PROJECT_DIR="/content/CFD_Solver_withCUDA"
COLAB_HOST="colab"  # SSH config alias (see setup instructions below)

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# =============================================================================
# Helper Functions
# =============================================================================

print_header() {
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}"
}

print_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

print_error() {
    echo -e "${RED}❌ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

print_info() {
    echo -e "${BLUE}ℹ️  $1${NC}"
}

# =============================================================================
# SSH Configuration Check
# =============================================================================

check_ssh_config() {
    SSH_CONFIG="$HOME/.ssh/config"
    
    if ! grep -q "Host $COLAB_HOST" "$SSH_CONFIG" 2>/dev/null; then
        print_warning "SSH config for Colab not found. Creating entry..."
        print_info "Please update the hostname after running Colab notebook"
        
        cat >> "$SSH_CONFIG" << EOF

# Google Colab SSH Configuration
Host $COLAB_HOST
    HostName your-cloudflared-hostname.trycloudflare.com
    User root
    Port 22
    StrictHostKeyChecking no
    UserKnownHostsFile /dev/null
    
EOF
        print_info "Created SSH config entry. Edit $SSH_CONFIG to update hostname"
        print_info "Replace 'your-cloudflared-hostname.trycloudflare.com' with actual hostname from Colab"
        return 1
    fi
    return 0
}

# =============================================================================
# Connection Test
# =============================================================================

test_connection() {
    print_info "Testing connection to Colab..."
    if ssh -o ConnectTimeout=5 "$COLAB_HOST" "echo 'Connection successful'" &>/dev/null; then
        print_success "Connected to Colab"
        return 0
    else
        print_error "Cannot connect to Colab"
        print_info "Make sure:"
        print_info "  1. Colab notebook is running"
        print_info "  2. SSH tunnel is active (Cell 2)"
        print_info "  3. SSH config hostname is correct"
        return 1
    fi
}

# =============================================================================
# Exclude Patterns
# =============================================================================

# Files/directories to exclude from sync
RSYNC_EXCLUDE=(
    "--exclude=build/"
    "--exclude=.git/"
    "--exclude=html/"
    "--exclude=latex/"
    "--exclude=.DS_Store"
    "--exclude=*.o"
    "--exclude=*.so"
    "--exclude=*.a"
    "--exclude=*.dylib"
    "--exclude=CMakeCache.txt"
    "--exclude=CMakeFiles/"
    "--exclude=cmake_install.cmake"
    "--exclude=Makefile"
    "--exclude=*.vtk"
    "--exclude=*.vtu"
    "--exclude=.vscode/"
    "--exclude=*.swp"
    "--exclude=*.swo"
    "--exclude=*~"
)

# =============================================================================
# Push: Local → Colab
# =============================================================================

push_to_colab() {
    print_header "Pushing Local Changes to Colab"
    
    if ! test_connection; then
        exit 1
    fi
    
    print_info "Syncing: $LOCAL_PROJECT_DIR → $REMOTE_PROJECT_DIR"
    
    # Perform rsync
    rsync -avz --progress \
        "${RSYNC_EXCLUDE[@]}" \
        "$LOCAL_PROJECT_DIR/" \
        "$COLAB_HOST:$REMOTE_PROJECT_DIR/"
    
    if [ $? -eq 0 ]; then
        print_success "Push complete!"
        
        # Show remote status
        print_info "Remote directory contents:"
        ssh "$COLAB_HOST" "ls -lh $REMOTE_PROJECT_DIR | head -20"
    else
        print_error "Push failed!"
        exit 1
    fi
}

# =============================================================================
# Pull: Colab → Local
# =============================================================================

pull_from_colab() {
    print_header "Pulling Changes from Colab to Local"
    
    if ! test_connection; then
        exit 1
    fi
    
    print_warning "This will overwrite local files with Colab versions!"
    read -p "Continue? (y/N) " -n 1 -r
    echo
    
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_info "Pull cancelled"
        exit 0
    fi
    
    print_info "Syncing: $REMOTE_PROJECT_DIR → $LOCAL_PROJECT_DIR"
    
    # Perform rsync
    rsync -avz --progress \
        "${RSYNC_EXCLUDE[@]}" \
        "$COLAB_HOST:$REMOTE_PROJECT_DIR/" \
        "$LOCAL_PROJECT_DIR/"
    
    if [ $? -eq 0 ]; then
        print_success "Pull complete!"
        print_info "Local files updated"
    else
        print_error "Pull failed!"
        exit 1
    fi
}

# =============================================================================
# Status Check
# =============================================================================

check_status() {
    print_header "Checking Sync Status"
    
    if ! test_connection; then
        exit 1
    fi
    
    print_info "Comparing local and remote directories..."
    
    # Check for differences
    rsync -avzn --delete \
        "${RSYNC_EXCLUDE[@]}" \
        "$LOCAL_PROJECT_DIR/" \
        "$COLAB_HOST:$REMOTE_PROJECT_DIR/" \
        | grep -E "^(sending|receiving|deleting)" > /tmp/sync_status.txt
    
    if [ -s /tmp/sync_status.txt ]; then
        print_warning "Differences found:"
        cat /tmp/sync_status.txt
        echo ""
        print_info "Run './sync_to_colab.sh push' to sync to Colab"
        print_info "Run './sync_to_colab.sh pull' to sync from Colab"
    else
        print_success "Directories are in sync!"
    fi
    
    rm -f /tmp/sync_status.txt
}

# =============================================================================
# Git Sync (Alternative Method)
# =============================================================================

git_sync() {
    print_header "Git-based Sync"
    print_info "Using git to sync changes..."
    
    # Local: commit and push
    cd "$LOCAL_PROJECT_DIR"
    
    if [ -n "$(git status --porcelain)" ]; then
        print_info "Local changes detected. Committing..."
        git add .
        git commit -m "Auto-sync from local at $(date '+%Y-%m-%d %H:%M:%S')"
        git push
        print_success "Local changes pushed to GitHub"
    else
        print_info "No local changes to commit"
    fi
    
    # Remote: pull
    if test_connection; then
        print_info "Pulling changes on Colab..."
        ssh "$COLAB_HOST" "cd $REMOTE_PROJECT_DIR && git pull"
        print_success "Colab updated from GitHub"
    fi
}

# =============================================================================
# Setup Instructions
# =============================================================================

show_setup() {
    print_header "Setup Instructions"
    
    cat << 'EOF'

To use this sync script, follow these steps:

1. RUN COLAB NOTEBOOK
   - Upload and run: Colab_VSCode_Setup.ipynb
   - Execute Cell 2 to get SSH command
   - Copy the hostname (e.g., random-name.trycloudflare.com)

2. UPDATE SSH CONFIG
   Edit: ~/.ssh/config
   
   Add or update:
   ```
   Host colab
       HostName YOUR-HOSTNAME.trycloudflare.com
       User root
       Port 22
       StrictHostKeyChecking no
       UserKnownHostsFile /dev/null
   ```

3. TEST CONNECTION
   ```
   ssh colab
   ```
   
   If successful, exit and proceed.

4. USE SYNC SCRIPT
   ```
   ./sync_to_colab.sh push    # Send local → Colab
   ./sync_to_colab.sh pull    # Get Colab → local
   ./sync_to_colab.sh status  # Check differences
   ./sync_to_colab.sh git     # Use git to sync
   ```

TIP: Update the hostname in ~/.ssh/config each time you start a new
     Colab session (hostname changes each time).

EOF
}

# =============================================================================
# Main Script
# =============================================================================

main() {
    # Check if rsync is available
    if ! command -v rsync &> /dev/null; then
        print_error "rsync not found. Please install it:"
        print_info "  Mac: brew install rsync"
        print_info "  Linux: sudo apt-get install rsync"
        exit 1
    fi
    
    # Parse command
    case "${1:-}" in
        push)
            push_to_colab
            ;;
        pull)
            pull_from_colab
            ;;
        status)
            check_status
            ;;
        git)
            git_sync
            ;;
        setup)
            show_setup
            ;;
        test)
            test_connection
            ;;
        *)
            echo "CFD Solver - Colab Sync Tool"
            echo ""
            echo "Usage: $0 {push|pull|status|git|setup|test}"
            echo ""
            echo "Commands:"
            echo "  push    - Upload local changes to Colab"
            echo "  pull    - Download Colab changes to local"
            echo "  status  - Check sync status"
            echo "  git     - Use git to sync (alternative method)"
            echo "  setup   - Show setup instructions"
            echo "  test    - Test connection to Colab"
            echo ""
            exit 1
            ;;
    esac
}

# Run main function
main "$@"
