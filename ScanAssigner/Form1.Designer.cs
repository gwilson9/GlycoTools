namespace ScanAssigner
{
    partial class Form1
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.label1 = new System.Windows.Forms.Label();
            this.uniprotGlycoDBFile = new System.Windows.Forms.TextBox();
            this.rawFilePathsBox = new System.Windows.Forms.Label();
            this.outputFolderPath = new System.Windows.Forms.TextBox();
            this.label3 = new System.Windows.Forms.Label();
            this.button1 = new System.Windows.Forms.Button();
            this.prgProgress = new System.Windows.Forms.ProgressBar();
            this.resultFilesBox = new System.Windows.Forms.ListBox();
            this.filePathsBox = new System.Windows.Forms.ListBox();
            this.label2 = new System.Windows.Forms.Label();
            this.SuspendLayout();
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(15, 18);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(107, 13);
            this.label1.TabIndex = 0;
            this.label1.Text = "Byonic Results Paths";
            // 
            // uniprotGlycoDBFile
            // 
            this.uniprotGlycoDBFile.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.uniprotGlycoDBFile.Location = new System.Drawing.Point(14, 274);
            this.uniprotGlycoDBFile.Name = "uniprotGlycoDBFile";
            this.uniprotGlycoDBFile.Size = new System.Drawing.Size(432, 20);
            this.uniprotGlycoDBFile.TabIndex = 3;
            // 
            // rawFilePathsBox
            // 
            this.rawFilePathsBox.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Right)));
            this.rawFilePathsBox.AutoSize = true;
            this.rawFilePathsBox.Location = new System.Drawing.Point(225, 18);
            this.rawFilePathsBox.Name = "rawFilePathsBox";
            this.rawFilePathsBox.Size = new System.Drawing.Size(78, 13);
            this.rawFilePathsBox.TabIndex = 2;
            this.rawFilePathsBox.Text = "Raw File Paths";
            // 
            // outputFolderPath
            // 
            this.outputFolderPath.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.outputFolderPath.Location = new System.Drawing.Point(14, 316);
            this.outputFolderPath.Name = "outputFolderPath";
            this.outputFolderPath.Size = new System.Drawing.Size(432, 20);
            this.outputFolderPath.TabIndex = 5;
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(11, 297);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(71, 13);
            this.label3.TabIndex = 4;
            this.label3.Text = "Output Folder";
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(172, 343);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(93, 21);
            this.button1.TabIndex = 6;
            this.button1.Text = "Assign Spectra";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click);
            // 
            // prgProgress
            // 
            this.prgProgress.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.prgProgress.Location = new System.Drawing.Point(18, 372);
            this.prgProgress.Name = "prgProgress";
            this.prgProgress.Size = new System.Drawing.Size(415, 22);
            this.prgProgress.TabIndex = 7;
            // 
            // resultFilesBox
            // 
            this.resultFilesBox.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.resultFilesBox.FormattingEnabled = true;
            this.resultFilesBox.HorizontalScrollbar = true;
            this.resultFilesBox.Location = new System.Drawing.Point(18, 34);
            this.resultFilesBox.Name = "resultFilesBox";
            this.resultFilesBox.Size = new System.Drawing.Size(201, 212);
            this.resultFilesBox.TabIndex = 8;
            // 
            // filePathsBox
            // 
            this.filePathsBox.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Right)));
            this.filePathsBox.FormattingEnabled = true;
            this.filePathsBox.HorizontalScrollbar = true;
            this.filePathsBox.Location = new System.Drawing.Point(229, 34);
            this.filePathsBox.Name = "filePathsBox";
            this.filePathsBox.Size = new System.Drawing.Size(206, 212);
            this.filePathsBox.TabIndex = 9;
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(12, 258);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(108, 13);
            this.label2.TabIndex = 10;
            this.label2.Text = "Uniprot Glyco DB File";
            // 
            // Form1
            // 
            this.AllowDrop = true;
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(450, 411);
            this.Controls.Add(this.label2);
            this.Controls.Add(this.filePathsBox);
            this.Controls.Add(this.resultFilesBox);
            this.Controls.Add(this.prgProgress);
            this.Controls.Add(this.button1);
            this.Controls.Add(this.outputFolderPath);
            this.Controls.Add(this.label3);
            this.Controls.Add(this.uniprotGlycoDBFile);
            this.Controls.Add(this.rawFilePathsBox);
            this.Controls.Add(this.label1);
            this.Name = "Form1";
            this.Text = "Scan Assigner";
            this.Load += new System.EventHandler(this.Form1_Load);
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.TextBox uniprotGlycoDBFile;
        private System.Windows.Forms.Label rawFilePathsBox;
        private System.Windows.Forms.TextBox outputFolderPath;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.Button button1;
        private System.Windows.Forms.ProgressBar prgProgress;
        private System.Windows.Forms.ListBox resultFilesBox;
        private System.Windows.Forms.ListBox filePathsBox;
        private System.Windows.Forms.Label label2;
    }
}

