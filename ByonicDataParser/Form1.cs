using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Data.SQLite;
using System.IO;
using System.Reflection;
using System.Threading;

namespace ByonicDataParser
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
            this.DragEnter += new DragEventHandler(Form1_DragEnter);
            this.DragDrop += new DragEventHandler(Form1_DragDrop);
        }

        private void byonicResultsFiles_SelectedIndexChanged(object sender, EventArgs e)
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

            if (files.Any(f => Path.GetExtension(f).Equals(".byrslt")))
                e.Effect = DragDropEffects.Link;
        }

        private void Form1_DragDrop(object sender, DragEventArgs e)
        {
            if (!e.Data.GetDataPresent(DataFormats.FileDrop))
                return;

            string[] files = e.Data.GetData(DataFormats.FileDrop) as string[];
            if (files == null)
                return;

            foreach (string file in files.Where(f => Path.GetExtension(f).Equals(".byrslt")))
            {
                byonicResultsFiles.Items.Add(file);

                if (String.IsNullOrEmpty(outputPath.Text))
                {
                    outputPath.Text = Path.GetDirectoryName(file);
                }       
            }
        }

        private void button1_Click(object sender, EventArgs e)
        {
            List<string> paths = new List<string>();

            foreach(string item in byonicResultsFiles.Items)
            {
                paths.Add(item);
            }

            var byonicDataParser = new byonicDataParser(paths, outputPath.Text);
            //link event handlers here
            byonicDataParser.highlightListItems = HandleHighlightListItems;
            byonicDataParser.finish = HandleFinished;

            Thread thread = new Thread(byonicDataParser.parse);            
            thread.Start();
            
        }

        private void HandleHighlightListItems(object sender, HighlightEventArgs e)
        {
            HighlightListItems(e.resultsFile, e.resultsFile);
        }

        private void HighlightListItems(string resultsFile, string rawFile)
        {
            if (InvokeRequired)
            {
                byonicResultsFiles.Invoke(new Action<string, string>(HighlightListItems), resultsFile, rawFile);
                return;
            }

            for (int i = 0; i < byonicResultsFiles.Items.Count; i++)
            {           
                if (resultsFile.Equals(byonicResultsFiles.Items[i]))
                {
                    byonicResultsFiles.SetSelected(i, true);
                }
            }
        }     

        private void HandleFinished(object sender, EventArgs e)
        {
            finished("");
        }

        private void finished(string name)
        {
            if (InvokeRequired)
            {
                byonicResultsFiles.Invoke(new Action<string>(finished), name);
                return;
            }

            for (int i = 0; i < byonicResultsFiles.Items.Count; i++)
            {
                byonicResultsFiles.SetSelected(i, false);
            }
        }
    }
}
