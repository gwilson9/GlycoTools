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
using System.Reflection;
using System.Threading;
using System.Diagnostics;

namespace ScanAssigner
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();

            this.DragEnter += new DragEventHandler(Form1_DragEnter);
            this.DragDrop += new DragEventHandler(Form1_DragDrop);
        }

        private void byonicResultsPath_TextChanged(object sender, EventArgs e)
        {

        }

        private void Form1_DragEnter(object sender, DragEventArgs e)
        {            
            e.Effect = DragDropEffects.None;

            if (!e.Data.GetDataPresent(DataFormats.FileDrop))
                return;

            string[] files = e.Data.GetData(DataFormats.FileDrop) as string[];
            if (files == null)
                return;

            if (files.Any(f => Path.GetExtension(f).Equals(".raw") || Path.GetExtension(f).Equals(".txt") || Path.GetExtension(f).Equals(".csv")))
                e.Effect = DragDropEffects.Link;
        }

        private void Form1_DragDrop(object sender, DragEventArgs e)
        {
            if (!e.Data.GetDataPresent(DataFormats.FileDrop))
                return;

            string[] files = e.Data.GetData(DataFormats.FileDrop) as string[];
            if (files == null)
                return;

            foreach (string file in files.Where(f => Path.GetExtension(f).Equals(".raw")))
            {
                //rawFilePath.Text = file;
                filePathsBox.Items.Add(file);
                
            }

            foreach (string file in files.Where(f => Path.GetExtension(f).Equals(".txt")))
            {
                //byonicResultsPath.Text = file;

                resultFilesBox.Items.Add(file);

                if (string.IsNullOrEmpty(outputFolderPath.Text))
                    outputFolderPath.Text = Path.GetDirectoryName(file);
            }

            foreach(string file in files.Where(f => Path.GetExtension(f).Equals(".csv")))
            {
                uniprotGlycoDBFile.Text = file;
            }
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            
        }

        private void button1_Click(object sender, EventArgs e)
        {           
            var scanAssigner = new ScanAssigner(filePathsBox.Items.OfType<string>().ToList(), resultFilesBox.Items.OfType<string>().ToList(), outputFolderPath.Text, uniprotGlycoDBFile.Text);
            scanAssigner.UpdateProgress += HandleUpdateProgress;
            prgProgress.Value = prgProgress.Minimum;
            scanAssigner.highlightListItems += HandleHighlightListItems;
            Task thread = new Task(scanAssigner.Start);
            thread.Start();                     
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

            prgProgress.Value = (int) (prgProgress.Maximum * progressValue);
        }

        private void HandleHighlightListItems(object sender, HighlightEventArgs e)
        {
            HighlightListItems(e.resultsFile, e.rawFile);
        }

        private void HighlightListItems( string resultsFile, string rawFile )
        {
            if (InvokeRequired)
            {
                resultFilesBox.Invoke(new Action<string, string>(HighlightListItems), resultsFile, rawFile);
                return;     
            }

            for(int i = 0; i < resultFilesBox.Items.Count; i++)
            {
                if (resultsFile.Equals(resultFilesBox.Items[i]))
                {
                    resultFilesBox.SetSelected(i, true);
                }
            }

            for(int i = 0; i < filePathsBox.Items.Count; i++)
            {
                if (rawFile.Equals(filePathsBox.Items[i]))
                {
                    filePathsBox.SetSelected(i, true);
                }
            }
        }
    }
}
