using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.IO;

namespace DeglycoDataBrowser
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();

            this.DragEnter += new DragEventHandler(Form1_DragEnter);
            this.DragDrop += new DragEventHandler(Form1_DragDrop);
        }


        private void Form1_DragEnter(object sender, DragEventArgs e)
        {
            e.Effect = DragDropEffects.None;

            if (!e.Data.GetDataPresent(DataFormats.FileDrop))
                return;

            string[] files = e.Data.GetData(DataFormats.FileDrop) as string[];
            if (files == null)
                return;

            if (files.Any(f => Path.GetExtension(f).Equals(".csv") || Path.GetExtension(f).Equals(".fasta")))
                e.Effect = DragDropEffects.Link;
        }

        private void Form1_DragDrop(object sender, DragEventArgs e)
        {
            if (!e.Data.GetDataPresent(DataFormats.FileDrop))
                return;

            string[] files = e.Data.GetData(DataFormats.FileDrop) as string[];
            if (files == null)
                return;

            foreach (string file in files.Where(f => Path.GetExtension(f).Equals(".csv")))
            {
                //rawFilePath.Text = file;
                peptideFiles.Items.Add(file);

                if (string.IsNullOrEmpty(outputFolder.Text))
                {
                    outputFolder.Text = Path.GetDirectoryName(file);
                }
            }

            foreach (string file in files.Where(f => Path.GetExtension(f).Equals(".fasta")))
            {
                uniprotDB.Text = file;               
            }
        }

        private void button1_Click(object sender, EventArgs e)
        {
            List<string> peptideFile = peptideFiles.Items.OfType<string>().ToList();
            string db = uniprotDB.Text;
            string outputPath = outputFolder.Text;

            deglycoDBGenerator dbGen = new deglycoDBGenerator(peptideFile, db, outputPath);
            dbGen.UpdateProgress += HandleUpdateProgress;
            dbGen.HighlightActiveFile += HandleHighlightListItems;
            Task task = new Task(dbGen.Start);
            task.Start(); 
        }

        private void peptideFiles_SelectedIndexChanged(object sender, EventArgs e)
        {

        }

        private void HandleUpdateProgress(object sender, ProgressEventArgs e)
        {
            ChangeProgressBarValue(e.Progress);
        }

        private void ChangeProgressBarValue(double progressValue)
        {
            if (InvokeRequired)
            {
                prgProgress.Invoke(new Action<double>(ChangeProgressBarValue), progressValue);
                return;
            }
            prgProgress.Value = (int)(prgProgress.Maximum * progressValue);
        }

        private void HandleHighlightListItems(object sender, HighlightEventArgs e)
        {
            HighlightListItems(e.resultsFile);
        }

        private void HighlightListItems(string resultsFile)
        {
            if (InvokeRequired)
            {
                peptideFiles.Invoke(new Action<string>(HighlightListItems), resultsFile);
                return;
            }

            for (int i = 0; i < peptideFiles.Items.Count; i++)
            {
                if (resultsFile.Equals(peptideFiles.Items[i]))
                {
                    peptideFiles.SetSelected(i, true);
                }
            }
        }
    }
}
