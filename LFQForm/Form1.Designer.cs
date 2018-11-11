namespace LFQForm
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
            this.byonicFiles = new System.Windows.Forms.ListBox();
            this.button1 = new System.Windows.Forms.Button();
            this.rawFiles = new System.Windows.Forms.ListBox();
            this.dataGridView1 = new System.Windows.Forms.DataGridView();
            this.outputPath = new System.Windows.Forms.TextBox();
            ((System.ComponentModel.ISupportInitialize)(this.dataGridView1)).BeginInit();
            this.SuspendLayout();
            // 
            // byonicFiles
            // 
            this.byonicFiles.AllowDrop = true;
            this.byonicFiles.FormattingEnabled = true;
            this.byonicFiles.ItemHeight = 24;
            this.byonicFiles.Location = new System.Drawing.Point(15, 24);
            this.byonicFiles.Margin = new System.Windows.Forms.Padding(6);
            this.byonicFiles.Name = "byonicFiles";
            this.byonicFiles.Size = new System.Drawing.Size(220, 268);
            this.byonicFiles.TabIndex = 0;
            this.byonicFiles.SelectedIndexChanged += new System.EventHandler(this.listBox1_SelectedIndexChanged);
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(162, 387);
            this.button1.Margin = new System.Windows.Forms.Padding(6);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(172, 50);
            this.button1.TabIndex = 1;
            this.button1.Text = "Start";
            this.button1.UseVisualStyleBackColor = true;
            // 
            // rawFiles
            // 
            this.rawFiles.AllowDrop = true;
            this.rawFiles.FormattingEnabled = true;
            this.rawFiles.ItemHeight = 24;
            this.rawFiles.Location = new System.Drawing.Point(266, 24);
            this.rawFiles.Name = "rawFiles";
            this.rawFiles.Size = new System.Drawing.Size(237, 268);
            this.rawFiles.TabIndex = 2;
            // 
            // dataGridView1
            // 
            this.dataGridView1.ColumnHeadersHeightSizeMode = System.Windows.Forms.DataGridViewColumnHeadersHeightSizeMode.AutoSize;
            this.dataGridView1.Location = new System.Drawing.Point(122, 483);
            this.dataGridView1.Name = "dataGridView1";
            this.dataGridView1.RowTemplate.Height = 31;
            this.dataGridView1.Size = new System.Drawing.Size(240, 150);
            this.dataGridView1.TabIndex = 3;
            // 
            // outputPath
            // 
            this.outputPath.Location = new System.Drawing.Point(15, 338);
            this.outputPath.Name = "outputPath";
            this.outputPath.Size = new System.Drawing.Size(490, 29);
            this.outputPath.TabIndex = 4;
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(11F, 24F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(528, 645);
            this.Controls.Add(this.outputPath);
            this.Controls.Add(this.dataGridView1);
            this.Controls.Add(this.rawFiles);
            this.Controls.Add(this.button1);
            this.Controls.Add(this.byonicFiles);
            this.Margin = new System.Windows.Forms.Padding(6);
            this.Name = "Form1";
            this.Text = "Form1";
            ((System.ComponentModel.ISupportInitialize)(this.dataGridView1)).EndInit();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.ListBox byonicFiles;
        private System.Windows.Forms.Button button1;
        private System.Windows.Forms.ListBox rawFiles;
        private System.Windows.Forms.DataGridView dataGridView1;
        private System.Windows.Forms.TextBox outputPath;
    }
}

