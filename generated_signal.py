#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 15:06:03 2024

@author: richy
"""
import numpy as np

# Function to generate and save a signal
def generate_signal(filepath):
    # Generate a simple sinusoidal signal
    t = np.linspace(0, 1, 1000)  # Time vector
    frequency = 5  # Frequency in Hz
    amplitude = 1  # Amplitude
    signal = amplitude * np.sin(2 * np.pi * frequency * t)

    # Save the signal to the specified file
    np.savetxt(filepath, signal)
    print(f"Signal saved to {filepath}")

if __name__ == "__main__":
    filepath = "generated_signal.txt"  # Specify the file path
    generate_signal(filepath)