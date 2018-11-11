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

            deglycoDBGenerator dbGen = new deglycoDBGenerator(peptideFile, db);
            Task task = new Task(dbGen.Start);
            task.Start();
        }
    }
}
