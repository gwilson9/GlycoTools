using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace GlycoTools
{
    public partial class GlycoToolsForm : Form
    {
        public GlycoToolsForm()
        {
            InitializeComponent();
        }

        private void toolStripButton1_Click(object sender, EventArgs e)
        {

        }

        private void tsbTool1_Click(object sender, EventArgs e)
        {
            var scanAssigner = new ScanAssigner.Form1{ MdiParent = this };
            scanAssigner.Show();
        }

        private void toolStripButton1_Click_1(object sender, EventArgs e)
        {
            var byonicResultsParser = new ByonicDataParser.Form1 { MdiParent = this };
            byonicResultsParser.Show();
        }

        private void toolStripButton2_Click(object sender, EventArgs e)
        {
            var glycoCompiler = new GlycoCompiler.Form1 { MdiParent = this };
            glycoCompiler.Show();
        }

        private void toolStripButton1_Click_2(object sender, EventArgs e)
        {
            var deglycoDBGenerator = new DeglycoDataBrowser.Form1 { MdiParent = this };
            deglycoDBGenerator.Show();
        }
    }
}
