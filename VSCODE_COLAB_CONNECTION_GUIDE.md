# VS Code to Google Colab Connection Guide

Complete guide for connecting VS Code to Google Colab for remote CUDA development.

---

## 📋 Prerequisites

### On Your Local Machine:
- ✅ Visual Studio Code installed
- ✅ Internet connection
- ✅ SSH client (built-in on Mac/Linux)

### On Google Colab:
- ✅ Free Google account
- ✅ GPU runtime enabled

---

## 🚀 Step-by-Step Setup

### **Part 1: Setup VS Code Extensions**

1. **Install Remote - SSH Extension**
   - Open VS Code
   - Press `Cmd+Shift+X` (Mac) or `Ctrl+Shift+X` (Windows/Linux)
   - Search: "Remote - SSH"
   - Install: `ms-vscode-remote.remote-ssh`

2. **Verify Installation**
   - Look for Remote icon in left sidebar (looks like `><`)
   - Or press `Cmd+Shift+P` and type "Remote-SSH" to see commands

---

### **Part 2: Setup Google Colab**

1. **Upload Notebook to Colab**
   - Go to https://colab.research.google.com
   - Click: `File` → `Upload notebook`
   - Upload: `Colab_VSCode_Setup.ipynb` (from this repository)

2. **Enable GPU**
   - Click: `Runtime` → `Change runtime type`
   - Hardware accelerator: Select `GPU`
   - GPU type: Choose `T4` (free tier) or better
   - Click `Save`

3. **Run Setup Cells**
   - Execute cells 1-5 in order
   - Wait for Cell 2 to display SSH command
   - **⚠️ IMPORTANT:** Copy the SSH command (looks like):
     ```
     ssh root@random-name.trycloudflare.com
     ```
   - Keep this tab open!

---

### **Part 3: Connect VS Code to Colab**

1. **Open Remote SSH**
   - In VS Code, press `Cmd+Shift+P` (Mac) or `Ctrl+Shift+P` (Windows)
   - Type: `Remote-SSH: Connect to Host`
   - Select it

2. **Enter SSH Command**
   - Paste the SSH command from Colab (Cell 2 output)
   - Press Enter

3. **Select Platform**
   - When prompted, select: `Linux`

4. **Enter Password**
   - Password: `cfd_solver_2026` (or what you set in the notebook)
   - Check "Remember password" for convenience

5. **Wait for Connection**
   - VS Code will install VS Code Server on Colab
   - This takes 1-2 minutes on first connection
   - Status shown in bottom-left corner

6. **Open Folder**
   - Click: `File` → `Open Folder`
   - Enter path: `/content/CFD_Solver_withCUDA`
   - Click `OK`

**🎉 You're now connected!** VS Code is running on Colab's GPU instance.

---

## 🔨 Building Your Project

### **Option 1: Using VS Code Terminal**

1. **Open Terminal**
   - In VS Code: `` Ctrl+` `` or `View` → `Terminal`

2. **Navigate and Build**
   ```bash
   cd build
   cmake .. -DCMAKE_BUILD_TYPE=Release
   make -j$(nproc)
   ```

3. **Run Your Solver**
   ```bash
   ./CFD_solver_gpu ../json_Files/Solver_Config.json
   ```

### **Option 2: Using VS Code Tasks**

Create `.vscode/tasks.json` in your project:

```json
{
  "version": "2.0.0",
  "tasks": [
    {
      "label": "Build CFD Solver GPU",
      "type": "shell",
      "command": "cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make -j$(nproc)",
      "group": {
        "kind": "build",
        "isDefault": true
      },
      "problemMatcher": ["$gcc"]
    },
    {
      "label": "Run CFD Solver GPU",
      "type": "shell",
      "command": "./build/CFD_solver_gpu ./json_Files/Solver_Config.json",
      "dependsOn": ["Build CFD Solver GPU"]
    }
  ]
}
```

Then press `Cmd+Shift+B` to build!

---

## 🐛 Troubleshooting

### **Connection Issues**

**Problem:** "Could not establish connection"
- **Solution:** 
  - Check Colab notebook is still running
  - Re-run Cell 2 in Colab to refresh SSH tunnel
  - Copy new SSH command and reconnect

**Problem:** "Host key verification failed"
- **Solution:**
  ```bash
  # On Mac/Linux terminal:
  ssh-keygen -R hostname.trycloudflare.com
  ```
  - Replace `hostname` with your actual hostname from SSH command

**Problem:** "Connection timeout"
- **Solution:**
  - Check your internet connection
  - Colab might be restarting; wait 1-2 minutes
  - Try re-running Colab setup cells

### **Build Issues**

**Problem:** "nvcc: command not found"
- **Solution:**
  ```bash
  # In VS Code terminal:
  export PATH=/usr/local/cuda/bin:$PATH
  ```

**Problem:** "No CMAKE_CUDA_COMPILER could be found"
- **Solution:**
  ```bash
  # Verify CUDA installation:
  which nvcc
  nvcc --version
  
  # If missing, re-run Colab Cell 1
  ```

**Problem:** "fatal error: json/json.h: No such file"
- **Solution:**
  ```bash
  # Install missing dependency:
  sudo apt-get install -y libjsoncpp-dev
  ```

### **Session Management**

**Problem:** Colab disconnects after inactivity
- **Solution:**
  - Run Cell 7 (Keep Alive) in Colab notebook
  - Keep Colab browser tab open
  - Colab has 12-hour max session limit (free tier)

**Problem:** Lost work after disconnection
- **Solution:**
  - Always commit changes to git:
    ```bash
    git add .
    git commit -m "Progress update"
    git push
    ```
  - Or sync to Google Drive (see next section)

---

## 💾 Syncing Files Between Local & Colab

### **Option 1: Git (Recommended)**

```bash
# In VS Code terminal (connected to Colab):
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"

# After making changes:
git add .
git commit -m "Your changes"
git push
```

### **Option 2: Automatic Sync Script**

Use the provided `sync_to_colab.sh` script:

```bash
# On your local machine:
./sync_to_colab.sh push   # Upload changes to Colab
./sync_to_colab.sh pull   # Download changes from Colab
```

### **Option 3: VS Code Sync Extension**

- Install: "Settings Sync" extension
- Sync your settings, extensions, and keybindings

---

## 📊 Monitoring GPU Usage

### **In VS Code Terminal:**

```bash
# Watch GPU utilization
watch -n 1 nvidia-smi

# Check CUDA version
nvcc --version

# Monitor during simulation
nvidia-smi --query-gpu=utilization.gpu,utilization.memory,temperature.gpu --format=csv -l 1
```

### **In Colab Notebook:**

Add a monitoring cell:
```python
!watch -n 1 nvidia-smi
```

---

## 🎯 Best Practices

1. **Regular Commits**
   - Commit every significant change
   - Colab sessions can disconnect unexpectedly

2. **Use Branches**
   ```bash
   git checkout -b colab-testing
   # Make experimental changes
   git push -u origin colab-testing
   ```

3. **Keep Notebooks Running**
   - Don't close Colab browser tab
   - Run keep-alive cell (Cell 7)

4. **Monitor Resource Usage**
   - Check GPU memory before large compilations
   - Free tier: ~12GB GPU RAM

5. **Save Outputs**
   - Copy VTK files to Google Drive
   - Download important results before session ends

---

## 🚦 Quick Reference Commands

### **Connection**
```bash
# From local VS Code
Cmd+Shift+P → "Remote-SSH: Connect to Host"
# Paste SSH command from Colab
```

### **Build**
```bash
cd /content/CFD_Solver_withCUDA/build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

### **Run**
```bash
./CFD_solver_gpu ../json_Files/Solver_Config.json
```

### **Disconnect**
```bash
# In VS Code
File → Close Remote Connection
# Or click on remote indicator (bottom-left)
```

---

## 📱 Mobile Development (Bonus)

You can even code on mobile using:
- **Code Server** on Colab
- Access via browser on phone/tablet

Setup in Colab:
```python
!curl -fsSL https://code-server.dev/install.sh | sh
!code-server --bind-addr 0.0.0.0:8080 --auth none &
!npm install -g localtunnel
!lt --port 8080
```

---

## 🆘 Need Help?

### **Resources:**
- VS Code Remote-SSH Docs: https://code.visualstudio.com/docs/remote/ssh
- Colab SSH GitHub: https://github.com/WassimBenzarti/colab-ssh
- CUDA Toolkit Docs: https://docs.nvidia.com/cuda/

### **Common Workflows:**

**Workflow 1: Quick Testing**
1. Connect to Colab via SSH
2. Edit code in VS Code
3. Build in integrated terminal
4. Run and observe results
5. Commit changes

**Workflow 2: Development**
1. Edit locally in VS Code (no connection)
2. Commit and push to GitHub
3. Connect to Colab
4. Pull latest changes
5. Build and test on GPU
6. Repeat

**Workflow 3: Continuous Development**
1. Keep VS Code connected to Colab
2. Edit files (auto-saved on Colab)
3. Build incrementally
4. Use git for checkpoints

---

## ⏰ Session Limits

### **Google Colab Free Tier:**
- **Maximum session:** 12 hours
- **Idle timeout:** 90 minutes (use keep-alive)
- **GPU allocation:** Not guaranteed during peak times

### **Google Colab Pro ($10/month):**
- **Maximum session:** 24 hours
- **Idle timeout:** Extended
- **GPU priority:** Higher allocation priority
- **Better GPUs:** V100, A100 access

---

## ✅ Connection Checklist

Before starting work:
- [ ] Colab notebook running with GPU enabled
- [ ] Cell 2 executed (SSH tunnel active)
- [ ] VS Code connected successfully
- [ ] Repository cloned at `/content/CFD_Solver_withCUDA`
- [ ] Build directory created
- [ ] CMakeLists configured for Colab
- [ ] CUDA verified (`nvcc --version`)
- [ ] Keep-alive activated (Cell 7)

---

## 📌 Pro Tips

1. **Use tmux on Colab**
   ```bash
   apt-get install tmux
   tmux new -s cfd
   # Your session persists if SSH disconnects
   ```

2. **Mount Google Drive for Persistent Storage**
   ```python
   from google.colab import drive
   drive.mount('/content/drive')
   # Save outputs to /content/drive/MyDrive/
   ```

3. **Speed Up Builds**
   ```bash
   # Use ccache
   apt-get install ccache
   export PATH=/usr/lib/ccache:$PATH
   ```

4. **Profile GPU Code**
   ```bash
   nvprof ./CFD_solver_gpu config.json
   # Or use Nsight Systems
   ```

---

**Happy Remote CUDA Development! 🚀**
