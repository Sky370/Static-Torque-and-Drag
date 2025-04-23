# Static Torque & Drag Model (v3.0)

A modular, 3D **soft-string static** torque and drag model for simulating drillstring behavior during tripping operations under varying circulation conditions. This model **excludes vibrational effects** and is tailored for practical field scenarios.

---

## 🚀 Features

- Simulate **Run-In-Hole (RIH)** and **Pull-Out-Of-Hole (POOH)** operations  
- Dedicated **torque analysis mode** ("TQ") for quick evaluations  
- Incorporates **viscous (hydraulic)** and **inertial fluid** forces  
- Supports **custom flow rate (GPM)** inputs that influence both hookload and torque  
- Includes **plotting utilities** for single and multi-rate comparisons  

---

## 🔧 What’s New in v3.0

- **Cleaned input workflow**: simpler arguments, explicit `GPM` control  
- **Modular plotting system**:
  - `lib/plots.py`: two-panel Hookload & Torque (single rate)
  - `lib/multi_plots.py`: compact multi-rate profile comparison  
- **Advanced pressure-loss modeling** in `lib/circulation.py`  
- **Hole vs. bit depth** toggle now supported (on/off-bottom logic)  
- **Refined structure**: 
  - Old scripts removed  
  - `.gitignore` added  
  - Unified entrypoint: `Run.py`  

> ℹ️ The fluid-force model is experimental. Validation with real field data is recommended.

---

## ⚙️ Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/Sky370/Torque-and-Drag.git
   cd Torque-and-Drag
2. **Create & activate a virtual environment**
    ```bash
    python3 -m venv .venv
    source .venv/bin/activate     # macOS/Linux
    .venv\Scripts\activate.bat    # Windows
3. **Install dependencies**
    ```bash
    pip install -r requirements.txt

---

## 📈 Usage

First, make sure the input data is imported or updated. A preloaded synthetic dataset (`NewData.xlsx`) is available in the `.Input/` directory.

Navigate to the `Run.py` file and adjust the arguments before executing the simulation.

### 🔹 Available Functions

#### `plot_hookload_torque(GPM=0)`
- Calculates hookload and torque **without fluid forces** (0 GPM)
- Runs calculations for all operations: `'POOH'`, `'RIH'`, and `'ROB'`

#### `plot_profiles_vs_gpm(GPM=[0, 250, 470], operation='ROB')`
- Calculates hookload and torque **with fluid circulation**
- Accepts multiple GPM values for comparison. Additional flowrates can be easily added/removed.
- Default operation is `'ROB'` (rotation off-bottom, no tripping)

---

## 🎛️ Argument Descriptions

- `bit_depth`: Bit depth in **feet**
- `wob`, `tob`: Weight on bit and torque on bit (**SI units**)
- `GPM`: Flow rate in **gallons per minute**
  - `plot_hookload_torque` takes a **single number**
  - `plot_profiles_vs_gpm` accepts a **list of numbers**
- `operation`: Choose from `['POOH', 'RIH', 'ROB']`
  - Not required for `plot_hookload_torque()` as it evaluates all operations

---

## 💾 Output Options

You can choose to save the simulation results in **interactive** or **static** formats:

- **Interactive (HTML)**  
  Uncomment the `html_out` argument to enable output in `.html` format

- **Static (PNG)**  
  Uncomment `image_out` and set `image_scale` for high-resolution `.png` images  
  - `image_scale` adjusts the output quality (higher value = higher resolution)

---
