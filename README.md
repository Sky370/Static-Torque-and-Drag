# Static Torque & Drag Model (v3.0)

A 3D soft-string **static** torque and drag model with no vibrational effects, suitable for simulating well tripping operations under a variety of circulation rates.

---

## Features

- **Tripping operations**: simulate Run-In-Hole (RIH) and Pull-Out-Of-Hole (POOH)  
- **Torque calculation** via dedicated “TQ” operation  
- **Viscous** (hydraulic) and **inertial** fluid forces fully accounted  
- **Configurable flow rate** (GPM) drives both hookload and torque profiles  
- **Modular plotting utilities** for single-rate and multi-rate comparison  

---

## What’s New in v3.0

- **Overhauled input workflow**: cleaner argument lists, explicit `GPM` parameter  
- **New plotting modules**  
  - `lib/plots.py` → single-rate, two-panel Hookload & Torque  
  - `lib/multi_plots.py` → multi-rate, single-panel “ROB vs. MD” (or any operation)  
- **Removed legacy files** and added a `.gitignore`  
- **Advanced pressure-loss** calculation in `lib/circulation.py`  
- **Hole-vs-bit depth toggle** (on/off-bottom) support  
- **New entrypoint**: `Run.py` replaces older scripts  
- **Fluid-force model**: please verify against field data  

---

## Installation

1. **Clone the repository**  
   ```bash
   git clone [https://github.com/Sky370/Torque-and-Drag.git](https://github.com/Sky370/Torque-and-Drag.git)
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
