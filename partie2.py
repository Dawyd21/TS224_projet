import tkinter as tk
from tkinter import filedialog, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class SignalProcessingApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Signal Processing Interface")

        # Buttons frame
        self.buttons_frame = tk.Frame(self.root)
        self.buttons_frame.pack(side=tk.LEFT, padx=10, pady=10)

        # Canvas frame
        self.canvas_frame = tk.Frame(self.root)
        self.canvas_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # Initialize signal variables
        self.signal = None
        self.bruited_signal = None

        # Buttons
        tk.Button(self.buttons_frame, text="Load Signal", command=self.load_signal).pack(pady=5)
        tk.Button(self.buttons_frame, text="Add Noise", command=self.add_noise).pack(pady=5)
        tk.Button(self.buttons_frame, text="Reload Signal", command=self.reload_signal).pack(pady=5)
        tk.Button(self.buttons_frame, text="Show Periodogram", command=self.show_periodogram).pack(pady=5)
        tk.Button(self.buttons_frame, text="Show Profile", command=self.show_profile).pack(pady=5)
        tk.Button(self.buttons_frame, text="Show Profile Parts", command=self.show_profile_parts).pack(pady=5)
        tk.Button(self.buttons_frame, text="Show Residual", command=self.show_residual).pack(pady=5)
        tk.Button(self.buttons_frame, text="Show F2(N)", command=self.show_f2n).pack(pady=5)

        # Placeholder for Matplotlib canvas
        self.canvas = None

    def load_signal(self):
        filepath = filedialog.askopenfilename(title="Select a Signal File", filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")])
        if filepath:
            try:
                self.signal = np.loadtxt(filepath)
                self.display_signal(self.signal, title="Loaded Signal")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load signal: {e}")

    def add_noise(self):
        if self.signal is None:
            messagebox.showwarning("Warning", "Please load a signal first.")
            return

        try:
            snr = float(filedialog.askstring("SNR", "Enter SNR (dB):"))
            power_signal = np.mean(self.signal ** 2)
            power_noise = power_signal / (10 ** (snr / 10))
            noise = np.sqrt(power_noise) * np.random.randn(len(self.signal))
            self.bruited_signal = self.signal + noise
            self.display_signal(self.bruited_signal, title=f"Bruited Signal (SNR: {snr} dB)")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to add noise: {e}")

    def reload_signal(self):
        if self.signal is None:
            messagebox.showwarning("Warning", "Please load a signal first.")
            return
        self.display_signal(self.signal, title="Reloaded Signal")

    def show_periodogram(self):
        if self.bruited_signal is None:
            messagebox.showwarning("Warning", "Please add noise to the signal first.")
            return
        freqs = np.fft.fftfreq(len(self.bruited_signal))
        spectrum = np.abs(np.fft.fft(self.bruited_signal)) ** 2
        self.display_signal(spectrum, x_axis=freqs, title="Periodogram")

    def show_profile(self):
        messagebox.showinfo("Info", "Functionality not yet implemented.")

    def show_profile_parts(self):
        messagebox.showinfo("Info", "Functionality not yet implemented.")

    def show_residual(self):
        messagebox.showinfo("Info", "Functionality not yet implemented.")

    def show_f2n(self):
        messagebox.showinfo("Info", "Functionality not yet implemented.")

    def display_signal(self, signal, x_axis=None, title="Signal"):
        if self.canvas:
            self.canvas.get_tk_widget().destroy()

        fig, ax = plt.subplots()
        if x_axis is None:
            x_axis = np.arange(len(signal))
        ax.plot(x_axis, signal)
        ax.set_title(title)

        self.canvas = FigureCanvasTkAgg(fig, master=self.canvas_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

if __name__ == "__main__":
    root = tk.Tk()
    app = SignalProcessingApp(root)
    root.mainloop()
