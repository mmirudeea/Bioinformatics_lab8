import tkinter as tk
from tkinter import ttk, messagebox

def align_dna(seq1, seq2, match_score, mismatch_score, gap_penalty):
    n, m = len(seq1), len(seq2)
    
    matrix = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    
    for i in range(n + 1):
        matrix[i][0] = i * gap_penalty
    for j in range(m + 1):
        matrix[0][j] = j * gap_penalty
        
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i - 1] == seq2[j - 1]:
                diag_score = matrix[i - 1][j - 1] + match_score
            else:
                diag_score = matrix[i - 1][j - 1] + mismatch_score
                
            up_score = matrix[i - 1][j] + gap_penalty
            left_score = matrix[i][j - 1] + gap_penalty
            
            matrix[i][j] = max(diag_score, up_score, left_score)
            
    align1, align2 = "", ""
    i, j = n, m
    
    while i > 0 and j > 0:
        current_score = matrix[i][j]
        diag = matrix[i - 1][j - 1]
        up = matrix[i - 1][j]
        
        score = match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score
            
        if current_score == diag + score:
            align1 += seq1[i - 1]
            align2 += seq2[j - 1]
            i -= 1
            j -= 1
        elif current_score == up + gap_penalty:
            align1 += seq1[i - 1]
            align2 += "-"
            i -= 1
        else:
            align1 += "-"
            align2 += seq2[j - 1]
            j -= 1
            
    while i > 0:
        align1 += seq1[i - 1]
        align2 += "-"
        i -= 1
    while j > 0:
        align1 += "-"
        align2 += seq2[j - 1]
        j -= 1
        
    return align1[::-1], align2[::-1], matrix[n][m]

class DNAAlignerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("DNA Sequence Aligner")
        self.root.geometry("600x500")
        
        input_frame = ttk.LabelFrame(root, text="Sequences")
        input_frame.pack(padx=10, pady=10, fill="x")
        
        ttk.Label(input_frame, text="Sequence 1:").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.seq1_entry = ttk.Entry(input_frame, width=50)
        self.seq1_entry.insert(0, "CGTGA")
        self.seq1_entry.grid(row=0, column=1, padx=5, pady=5)
        
        ttk.Label(input_frame, text="Sequence 2:").grid(row=1, column=0, padx=5, pady=5, sticky="w")
        self.seq2_entry = ttk.Entry(input_frame, width=50)
        self.seq2_entry.insert(0, "ACGTC")
        self.seq2_entry.grid(row=1, column=1, padx=5, pady=5)
        
        param_frame = ttk.LabelFrame(root, text="Parameters")
        param_frame.pack(padx=10, pady=10, fill="x")
        
        ttk.Label(param_frame, text="Match Score:").grid(row=0, column=0, padx=5, pady=5)
        self.match_var = tk.IntVar(value=1)
        ttk.Entry(param_frame, textvariable=self.match_var, width=10).grid(row=0, column=1, padx=5, pady=5)
        
        ttk.Label(param_frame, text="Mismatch Penalty:").grid(row=0, column=2, padx=5, pady=5)
        self.mismatch_var = tk.IntVar(value=-1)
        ttk.Entry(param_frame, textvariable=self.mismatch_var, width=10).grid(row=0, column=3, padx=5, pady=5)
        
        ttk.Label(param_frame, text="Gap Penalty:").grid(row=0, column=4, padx=5, pady=5)
        self.gap_var = tk.IntVar(value=-1)
        ttk.Entry(param_frame, textvariable=self.gap_var, width=10).grid(row=0, column=5, padx=5, pady=5)
        
        align_btn = ttk.Button(root, text="Align Sequences", command=self.run_alignment)
        align_btn.pack(pady=10)
        
        result_frame = ttk.LabelFrame(root, text="Alignment Result")
        result_frame.pack(padx=10, pady=10, fill="both", expand=True)
        
        self.result_text = tk.Text(result_frame, font=("Courier", 12), height=10, bg="#f4f4f4")
        self.result_text.pack(padx=5, pady=5, fill="both", expand=True)
        self.result_text.config(state=tk.DISABLED)

    def run_alignment(self):
        s1 = self.seq1_entry.get().strip().upper()
        s2 = self.seq2_entry.get().strip().upper()
        
        if not s1 or not s2:
            messagebox.showerror("Error", "Please enter both sequences.")
            return
            
        try:
            match = self.match_var.get()
            mismatch = self.mismatch_var.get()
            gap = self.gap_var.get()
        except ValueError:
            messagebox.showerror("Error", "Scores must be integers.")
            return

        aligned_1, aligned_2, final_score = align_dna(s1, s2, match, mismatch, gap)
        
        match_line = "".join(["|" if a == b else " " for a, b in zip(aligned_1, aligned_2)])
        
        self.result_text.config(state=tk.NORMAL)
        self.result_text.delete(1.0, tk.END)
        self.result_text.insert(tk.END, f"Score: {final_score}\n\n")
        self.result_text.insert(tk.END, f"Seq 1: {aligned_1}\n")
        self.result_text.insert(tk.END, f"       {match_line}\n")
        self.result_text.insert(tk.END, f"Seq 2: {aligned_2}\n")
        self.result_text.config(state=tk.DISABLED)

if __name__ == "__main__":
    root = tk.Tk()
    app = DNAAlignerApp(root)
    root.mainloop()