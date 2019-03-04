namespace DeglycoDataBrowser
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
            this.uniprotDB = new System.Windows.Forms.TextBox();
            this.label1 = new System.Windows.Forms.Label();
            this.label2 = new System.Windows.Forms.Label();
            this.peptideFiles = new System.Windows.Forms.ListBox();
            this.button1 = new System.Windows.Forms.Button();
            this.label3 = new System.Windows.Forms.Label();
            this.outputFolder = new System.Windows.Forms.TextBox();
            this.prgProgress = new System.Windows.Forms.ProgressBar();
            this.SuspendLayout();
            // 
            // uniprotDB
            // 
            this.uniprotDB.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.uniprotDB.Location = new System.Drawing.Point(15, 226);
            this.uniprotDB.Name = "uniprotDB";
            this.uniprotDB.Size = new System.Drawing.Size(338, 20);
            this.uniprotDB.TabIndex = 0;
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(12, 208);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(90, 13);
            this.label1.TabIndex = 1;
            this.label1.Text = "Uniprot Database";
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(12, 9);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(210, 13);
            this.label2.TabIndex = 2;
            this.label2.Text = "Compass Peptides File ( from ProteinHoard)";
            // 
            // peptideFiles
            // 
            this.peptideFiles.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.peptideFiles.FormattingEnabled = true;
            this.peptideFiles.Location = new System.Drawing.Point(15, 28);
            this.peptideFiles.Name = "peptideFiles";
            this.peptideFiles.Size = new System.Drawing.Size(338, 173);
            this.peptideFiles.TabIndex = 3;
            this.peptideFiles.SelectedIndexChanged += new System.EventHandler(this.peptideFiles_SelectedIndexChanged);
            // 
            // button1
            // 
            this.button1.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.button1.Location = new System.Drawing.Point(92, 296);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(164, 27);
            this.button1.TabIndex = 4;
            this.button1.Text = "Start";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click);
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(12, 249);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(71, 13);
            this.label3.TabIndex = 6;
            this.label3.Text = "Output Folder";
            // 
            // outputFolder
            // 
            this.outputFolder.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.outputFolder.Location = new System.Drawing.Point(15, 267);
            this.outputFolder.Name = "outputFolder";
            this.outputFolder.Size = new System.Drawing.Size(338, 20);
            this.outputFolder.TabIndex = 5;
            // 
            // prgProgress
            // 
            this.prgProgress.Location = new System.Drawing.Point(15, 333);
            this.prgProgress.Name = "prgProgress";
            this.prgProgress.Size = new System.Drawing.Size(335, 23);
            this.prgProgress.TabIndex = 7;
            // 
            // Form1
            // 
            this.AllowDrop = true;
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(362, 368);
            this.Controls.Add(this.prgProgress);
            this.Controls.Add(this.label3);
            this.Controls.Add(this.outputFolder);
            this.Controls.Add(this.button1);
            this.Controls.Add(this.peptideFiles);
            this.Controls.Add(this.label2);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.uniprotDB);
            this.FormBorderStyle = System.Windows.Forms.FormBorderStyle.FixedSingle;
            this.Name = "Form1";
            this.Text = "Database Generator";
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.TextBox uniprotDB;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.ListBox peptideFiles;
        private System.Windows.Forms.Button button1;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.TextBox outputFolder;
        private System.Windows.Forms.ProgressBar prgProgress;
    }
}

