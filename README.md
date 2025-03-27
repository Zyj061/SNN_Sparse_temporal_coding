# 🧠 Brain-Inspired Perception Information Processing System

**Author:** Yajing Zheng  
**Date:** April 2017  

> 📄 This repository contains the **implementation code** for our paper:  
> **"Spike-Based Motion Estimation for Object Tracking Through Bio-Inspired Unsupervised Learning"**  
> 📰 Published in *IEEE Transactions on Circuits and Systems I: Regular Papers*  
> 🔗 [Read the paper on IEEE Xplore](https://ieeexplore.ieee.org/abstract/document/9985998)

---

## ⚠️ Quick Note

Since **Step 1** is time- and memory-consuming, we have saved the execution results in `DATA.mat`.  
The **true innovation** of this project lies in **Step 2**, so you can directly load the `DATA.mat` file in MATLAB to proceed with the computations.

---

## 🔬 1. Unsupervised Training with STDP-based HMAX

1. Run the script `stdpRum.m` located in the `STDP/script` directory to perform unsupervised training on the 3D dataset and extract image features.

   After the training, save all feature data into `DATA.mat`:

   ```matlab
   save 'DATA.mat'
   ```

---

## 🧠 2. Encoding and Supervised Training

1. Load the pre-saved data into the workspace:

   ```matlab
   load 'DATA.mat'
   ```

2. Run `OBJREGmain.m` for training and evaluation:

   - **`arithEnc`**: Arithmetic-based encoding 🧮  
     Inputs: `COMMON.firingTime`, `COMMON.localFiringSpike`, `COMMON.firingSpike`, and `COMMON.picScale`.  
     Output: `PtnTrSet` (spatiotemporal patterns for training) and `TrainLabels`.

   - **`TrainSnglN`**: Trains a weight matrix for each object. 💪  
     After training, the weight matrix is stored.

   - **`Tesing`**: Uses trained weights to classify test images. 🧪  
     Outputs:  
     `TeFiring` (output spike train),  
     `DistTePtns` (vR distance to target spike train).

   - **`ResultAnalysis`**: Evaluates recognition results based on `DistTePtns`. 📊

---

## 🙌 Support & Citation

If you find this project helpful, please consider giving it a ⭐️ or **citing our paper** to support our research!  
Your support encourages us to continue exploring brain-inspired AI. 🌟🧠

### 📚 Cite our work:

```bibtex
@article{zheng2022spike,
  title={Spike-Based Motion Estimation for Object Tracking Through Bio-Inspired Unsupervised Learning},
  author={Zheng, Yajing and others},
  journal={IEEE Transactions on Circuits and Systems I: Regular Papers},
  year={2022},
  publisher={IEEE}
}
```

🔗 [https://ieeexplore.ieee.org/abstract/document/9985998](https://ieeexplore.ieee.org/abstract/document/9985998)
