# 🔌 Connect VS Code to Google Colab via SSH

## ✅ Configuration Complete!

Your SSH connection has been configured with:
- **Host:** focused-left-plugin-guards.trycloudflare.com
- **Password:** cfd_solver_2026
- **SSH Alias:** colab-cfd

---

## 🚀 Connect Now (3 Steps):

### **Method 1: Quick Connect (Recommended)**

1. **Press:** `Cmd+Shift+P` (Mac) or `Ctrl+Shift+P` (Windows/Linux)
2. **Type:** `Remote-SSH: Connect to Host`
3. **Select:** `colab-cfd` from the list
4. **Enter password:** `cfd_solver_2026`

**Done!** VS Code will connect to Colab's GPU instance.

---

### **Method 2: Manual SSH Command**

If you prefer command line:

```bash
ssh root@focused-left-plugin-guards.trycloudflare.com
# Password: cfd_solver_2026
```

---

## 📁 After Connecting:

Once connected, you'll see "SSH: colab-cfd" in the bottom-left corner.

1. **Open Folder:**
   - Click: `File` → `Open Folder`
   - Enter: `/content/CFD_Solver_withCUDA`
   - Click: OK

2. **Open Terminal:**
   - Press: `` Ctrl+` `` or `View` → `Terminal`
   - You'll have direct access to Colab's GPU environment

3. **Build Your Project:**
   ```bash
   cd /content/CFD_Solver_withCUDA
   mkdir -p build && cd build
   cmake .. -DCMAKE_BUILD_TYPE=Release
   make -j$(nproc)
   ```

4. **Verify GPU:**
   ```bash
   nvidia-smi
   nvcc --version
   ```

---

## 🔧 Troubleshooting:

### **Issue: "Could not establish connection"**
- Check that Colab notebook is still running
- Verify the cloudflare tunnel is active
- Try reconnecting

### **Issue: "Host key verification failed"**
```bash
ssh-keygen -R focused-left-plugin-guards.trycloudflare.com
```

### **Issue: "Connection timeout"**
- The tunnel URL changes each Colab session
- Update the hostname in `~/.ssh/config`:
  ```bash
  nano ~/.ssh/config
  # Update HostName to new URL
  ```

---

## ⏰ Session Notes:

- **Cloudflare URL changes** each time you restart Colab
- **Session limit:** 12 hours (free tier)
- **Keep Colab tab open** while using SSH
- Update SSH config when URL changes

---

## 🎯 Quick Reference:

**Connect:**
```
Cmd+Shift+P → Remote-SSH: Connect to Host → colab-cfd
Password: cfd_solver_2026
```

**Disconnect:**
```
File → Close Remote Connection
```

**Reconnect:**
```
Click on "SSH: colab-cfd" in bottom-left → Connect to Host
```

---

## ✅ Success Indicators:

After connecting, you should see:
- Bottom-left shows: `SSH: colab-cfd`
- Terminal shows: `root@...` prompt
- `nvidia-smi` shows Tesla T4/V100/A100
- `/content/` directory accessible

---

## 🚀 Ready to Code!

You now have full VS Code access to Colab's GPU environment:
- ✅ Full terminal access
- ✅ IntelliSense and debugging
- ✅ GPU acceleration
- ✅ File editing with VS Code features
- ✅ Git integration

**Start connecting now using the steps above!** 🎉
