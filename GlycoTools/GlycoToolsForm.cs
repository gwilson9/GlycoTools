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

        private void GlycoToolsForm_Load(object sender, EventArgs e)
        {

        }
        // Database generator botton
        private void button2_Click(object sender, EventArgs e)
        {
            //this.BackColor = Color.FromArgb(0, 0, 0);
            var deglycoDBGenerator = new DeglycoDataBrowser.Form1 { MdiParent = this };

            deglycoDBGenerator.Show();

            panel1.Hide();
            richTextBox1.Text = "Database searching of intact glycopeptide analyses are benefited by" +
                                " using an informed search space (Khatri et al. (2017) Anal. Bioanal. " +
                                "Chem.). We use a filtered database consisting of proteins with validated" +
                                " deglycosylated sites (N à D at NxS/T) to search intact glycopeptides. " +
                                "‘Database Generator’ is designed to take as input the ‘parsimony_peptides’ " +
                                "file from a ‘Protein Hoard’ performed within the Coon OMSSA Proteomic Software " +
                                "Suite (COMPASS, Wenger et al. (2011) Proteomics) and the corresponding Uniprot " +
                                "FASTA file. The program will generate as output an annotated peptides file and " +
                                "a FASTA file consisting of only proteins for which sites of deglycosylation were " +
                                "found. An example file can be found within '..GlycoTools/Example Data/Database " +
                                "Generator/parsimony_peptides_filtered.csv'.";
            richTextBox1.Show();            
        }

        private void panel1_Paint(object sender, PaintEventArgs e)
        {

        }

        // Byonic Parser
        private void button3_Click(object sender, EventArgs e)
        {
            var byonicResultsParser = new ByonicDataParser.Form1 { MdiParent = this };
            byonicResultsParser.Show();
            panel1.Hide();
            richTextBox1.Text = "Glycopeptides from Byonic searches are contained within a '.byrslt' file. "+
                                "The Byonic Extract program will parse the file and export relevant information "+
                                "for each PSM to be used in subsequent processing methods.";
            richTextBox1.Show();
            
        }

        // Scan Assigner
        private void button1_Click(object sender, EventArgs e)
        {
            var scanAssigner = new ScanAssigner.Form1 { MdiParent = this };
            scanAssigner.Show();
            panel1.Controls.Remove(richTextBox1);
            panel1.Hide();
            richTextBox1.Text = "'Scan Profiler' will provide annotation and quantitation for the results from individual " +
                                "raw files. We search for glycan fragments, glycopeptide fragments, and oxonium ions. "+
                                "Intensity and m/z are provided for each observed ion to assist in manual spectrum annotation. " +
                                "We also draw from Uniprot site annotation files to provide information for known glycosites. "+
                                "Quantitation is performed during this step and is reported as the area of "+
                                "the eluting peak.";
            richTextBox1.Show();
            
        }

        // Glyco Compiler
        private void button4_Click(object sender, EventArgs e)
        {
            var glycoCompiler = new GlycoCompiler.Form1 { MdiParent = this };
            glycoCompiler.Show();
            panel1.Hide();
            richTextBox1.Text = "The compiler tool will aggregate data from input list(s) of PSMs into organized files of " +
                                "glycopeptides, glycoproteins, glycosites, and glycans. This function is particularly useful for "+
                                "aggregating data from across injection replicates.";
            richTextBox1.Show();
            
        }

        // Glyco Tool title
        private void label2_Click(object sender, EventArgs e)
        {
            panel1.Show();
            richTextBox1.Hide();
            
        }

        private void richTextBox1_TextChanged(object sender, EventArgs e)
        {

        }

        // Database generator side button
        private void button5_Click(object sender, EventArgs e)
        {
            var deglycoDBGenerator = new DeglycoDataBrowser.Form1 { MdiParent = this };
            deglycoDBGenerator.Show();
            panel1.Hide();
            richTextBox1.Text = "Database searching of intact glycopeptide analyses are benefited by" +
                                " using an informed search space (Khatri et al. (2017) Anal. Bioanal. " +
                                "Chem.). We use a filtered database consisting of proteins with validated" +
                                " deglycosylated sites (N à D at NxS/T) to search intact glycopeptides. " +
                                "‘Database Generator’ is designed to take as input the ‘parsimony_peptides’ " +
                                "file from a ‘Protein Hoard’ performed within the Coon OMSSA Proteomic Software " +
                                "Suite (COMPASS, Wenger et al. (2011) Proteomics) and the corresponding Uniprot " +
                                "FASTA file. The program will generate as output an annotated peptides file and " +
                                "a FASTA file consisting of only proteins for which sites of deglycosylation were " +
                                "found. An example file can be found within '..GlycoTools/Example Data/Database " +
                                "Generator/parsimony_peptides_filtered.csv'.";
            richTextBox1.Show();
        }

        // Byonic parser side button
        private void button6_Click(object sender, EventArgs e)
        {
            var byonicResultsParser = new ByonicDataParser.Form1 { MdiParent = this };
            byonicResultsParser.Show();
            panel1.Hide();
            richTextBox1.Text = "Glycopeptides from Byonic searches are contained within a '.byrslt' file. " +
                                "The Byonic Extract program will parse the file and export relevant information " +
                                "for each PSM to be used in subsequent processing methods.";
            richTextBox1.Show();
        }

        // Scan assigner side button
        private void button7_Click(object sender, EventArgs e)
        {
            var scanAssigner = new ScanAssigner.Form1 { MdiParent = this };
            scanAssigner.Show();
            panel1.Hide();
            richTextBox1.Text = "'Scan Profiler' will provide annotation and quantitation for the results from individual " +
                                "raw files. We search for glycan fragments, glycopeptide fragments, and oxonium ions. " +
                                "Intensity and m/z are provided for each observed ion to assist in manual spectrum annotation. " +
                                "We also draw from Uniprot site annotation files to provide information for known glycosites. " +
                                "Quantitation is performed during this step and is reported as the area of " +
                                "the eluting peak area.";
            richTextBox1.Show();
            
        }

        // Glyco compiler side button
        private void button8_Click(object sender, EventArgs e)
        {
            var glycoCompiler = new GlycoCompiler.Form1 { MdiParent = this };
            glycoCompiler.Show();
            panel1.Hide();
            richTextBox1.Text = "The compiler tool will aggregate data from input list(s) of PSMs into organized files of " +
                                "glycopeptides, glycoproteins, glycosites, and glycans. This function is particularly useful for " +
                                "aggregating data from across injection replicates.";
            richTextBox1.Show();
            
        }

        private void panel2_Paint(object sender, PaintEventArgs e)
        {

        }

        private void button5_Click_1(object sender, EventArgs e)
        {
            System.Diagnostics.Process.Start("http://coonlabs.com/");
        }
    }
}
