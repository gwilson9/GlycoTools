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

namespace GlycoCompiler
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

            if (files.Any(f => Path.GetExtension(f).Equals(".txt") || Path.GetExtension(f).Equals(".csv")))
                e.Effect = DragDropEffects.Link;
        }

        private void Form1_DragDrop(object sender, DragEventArgs e)
        {
            if (!e.Data.GetDataPresent(DataFormats.FileDrop))
                return;

            string[] files = e.Data.GetData(DataFormats.FileDrop) as string[];
            if (files == null)
                return;

            foreach (string file in files.Where(f => Path.GetExtension(f).Equals(".txt")))
            {
                listBox1.Items.Add(file);

                if (string.IsNullOrEmpty(textBox2.Text))
                    textBox2.Text = Path.GetDirectoryName(file);
            }

            foreach (string file in files.Where(f => Path.GetExtension(f).Equals(".csv")))
            {
                textBox1.Text = file;

                if (string.IsNullOrEmpty(textBox2.Text))
                    textBox2.Text = Path.GetDirectoryName(file);
            }
        }

        private void label1_Click(object sender, EventArgs e)
        {

        }

        private void listBox1_SelectedIndexChanged(object sender, EventArgs e)
        {

        }

        private void button1_Click(object sender, EventArgs e)
        {
            List<string> results = new List<string>();

            foreach(string file in listBox1.Items)
            {
                results.Add(file);
            }

            var glycoCompiler = new glycoCompiler(results, textBox1.Text, textBox2.Text);
            glycoCompiler.UpdateProgress += HandleUpdateProgress;
            prgProgress.Value = prgProgress.Minimum;

            Thread thread = new Thread(glycoCompiler.compile);
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
            prgProgress.Value = (int)(prgProgress.Maximum * progressValue);
        }
    }
}
