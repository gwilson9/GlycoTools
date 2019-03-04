using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using LumenWorks.Framework.IO.Csv;
using CSMSL;
using CSMSL.IO;
using CSMSL.Chemistry;
using CSMSL.Proteomics;
using System.Text.RegularExpressions;

namespace DeglycoDataBrowser
{
    class deglycoDBGenerator
    {
        public List<string> files;
        public string database;
        public string outputPath;

        public deglycoDBGenerator(List<string> files, string databasePath, string outputPath)
        {
            this.files = files;
            this.database = databasePath;
            this.outputPath = outputPath;

        }

        public void Start()
        {
            Dictionary<string, string> db = parseDB(database);

            foreach (string file in files)
            {
                OnHighlightActiveFile(file);

                OnUpdateProgress(0.5);

                File parsedFile = new File(file, db, Path.GetFileNameWithoutExtension(file));

                var resultsPath = outputPath + "\\" + Path.GetFileNameWithoutExtension(file) + "_parsedByGlycoTools.csv";

                printData(parsedFile, resultsPath, db);

                OnUpdateProgress(0.75);

                var fastaPath = outputPath + "\\" + Path.GetFileNameWithoutExtension(file) + "_LocalizedPSMs.fasta";

                printFastaFile(fastaPath, parsedFile.localizedDeglycoPSMsWithMotif, db);

                OnUpdateProgress(1.0);

                OnUpdateProgress(0.0);
            }            
        }

        public static void printData(File file, string pathout, Dictionary<string, string> db)
        {           

            StreamWriter writer = new StreamWriter(pathout);

            string firstLine = "Peptide,Protein,Mass,Charge,RawFile,Spectrum Number,Is Deamidated PSM,Mod Position Within Protein,Deglyco Mod Count,Sequon Amino Acid (S/T),Is Localized";
            writer.WriteLine(firstLine);

            foreach(psms psm in file.allPsms)
            {
                string nextLine = psm.sequence + "," + psm.protein + "," + psm.mass + "," + psm.charge + "," + psm.rawFileName + "," + psm.spectrumNumber + ",";

                if(psm.deglycoMods.Count() > 0)
                {
                    nextLine += "TRUE,";
                }
                else
                {
                    nextLine += "FALSE,";
                }

                for( int i = 0; i < psm.deglycoMods.Count(); i++ )
                {
                    nextLine += psm.deglycoMods[i];
                    if(i != psm.deglycoMods.Count - 1)
                    {
                        nextLine += ";";
                    }
                }

                nextLine += "," + psm.deglycoMods.Count() + ",";

                for( int i = 0; i < psm.deglycoMods.Count(); i++)
                {
                    // Position +2 from modifed asparagine

                    if(psm.deglycoMods[i] + 1 < db[">" + psm.protein].Length)
                        nextLine += db[">" + psm.protein][psm.deglycoMods[i]+1];

                    if(i != psm.deglycoMods.Count() - 1)
                    {
                        nextLine += ";";
                    }
                }

                nextLine += ",";

                if(file.localizedDeglycoPSMsWithMotif.Contains(psm, new PsmComparer()))
                {
                    nextLine += "TRUE";
                }
                else
                {
                    nextLine += "FALSE";
                }

                writer.WriteLine(nextLine);
            }
            writer.Close();
        }

        public static void printFastaFile(string pathout, List<psms> psms, Dictionary<string, string> db)
        {

            HashSet<string> fastaHeaders = new HashSet<string>();

            foreach (psms psm in psms)
            {
                fastaHeaders.Add(psm.protein);
            }

            StreamWriter writer = new StreamWriter(pathout);

            foreach (string id in fastaHeaders)
            {
                writer.WriteLine(">" + id);

                char[] sequence = db[">" + id].ToCharArray();

                string sequenceChunk = "";

                for (int i = 0; i < sequence.Length; i++)
                {
                    sequenceChunk += sequence[i];
                    if ((i + 1) % 80 == 0)
                    {
                        writer.WriteLine(sequenceChunk);
                        sequenceChunk = "";
                    }
                }
                if (!sequenceChunk.Equals(""))
                {
                    writer.WriteLine(sequenceChunk);
                }
            }

            writer.Close();

        }

        public static List<string> getSites(List<string> data, Dictionary<string, string> db)
        {
            Dictionary<string, string> uniprotIDs = new Dictionary<string, string>();

            foreach (KeyValuePair<string, string> pair in db)
            {
                string uniprotID = pair.Key.Split('|')[1];
                uniprotIDs.Add(uniprotID, pair.Key);
            }


            List<string> returnList = new List<string>();

            foreach (string line in data)
            {

                string line2 = line + ",";
                // line = Protein (fasta header), modified seqeunce, 
                string[] lineArray = line.Split(',');
                string prot = lineArray[0];

                string peptideMod = lineArray[1].Replace("N[+1]", "n");
                peptideMod = peptideMod.Replace("C[+57]", "C");
                peptideMod = peptideMod.Replace("Q[-17]", "Q");
                peptideMod = peptideMod.Replace("M[+16]", "M");
                peptideMod = peptideMod.Replace("F[+42]", "F");

                if (db.ContainsKey(uniprotIDs[prot]))
                {
                    if (Regex.Matches(db[uniprotIDs[prot]], peptideMod.ToUpper()).Count > 0)
                    {
                        foreach (Match m in Regex.Matches(db[uniprotIDs[prot]], peptideMod.ToUpper()))
                        {
                            int startPosition = m.Index;
                            foreach (Match n in Regex.Matches(peptideMod, "n"))
                            {
                                int site = startPosition + n.Index + 1;

                                returnList.Add(line2 + "N(" + site + ")");
                            }
                        }
                    }
                    else
                    {
                        Console.WriteLine("Could not match peptide " + peptideMod.ToUpper() + " with protein " + prot);
                    }

                    //returnList.Add(line2.TrimEnd(';'));                   
                }
                else
                {
                    Console.WriteLine("db does not contain " + prot);
                }
            }
            return returnList;
        }

        public static Dictionary<string, string> parseDB(string pathIn)
        {
            Dictionary<string, string> returnDict = new Dictionary<string, string>();

            List<string> badsequences = new List<string>();

            //StreamReader reader = new StreamReader(pathIn);

            var reader = new FastaReader(pathIn);

            foreach (var fasta in reader.ReadNextFasta())
            {
                if (!fasta.IsDecoy)
                {
                    string desc = fasta.Description.Replace(",", ";");
                    returnDict.Add(">" + desc, fasta.Sequence);
                }                    
            }
            return returnDict;
        }

        public event EventHandler<ProgressEventArgs> UpdateProgress;

        protected virtual void OnUpdateProgress(double progress)
        {
            var handler = UpdateProgress;
            if (handler != null)
            {
                handler(this, new ProgressEventArgs(progress));
            }
        }

        public event EventHandler<HighlightEventArgs> HighlightActiveFile;

        protected virtual void OnHighlightActiveFile(string file)
        {
            var handler = HighlightActiveFile;
            if(handler != null)
            {
                handler(this, new HighlightEventArgs(file));
            }

        }

    }
}
