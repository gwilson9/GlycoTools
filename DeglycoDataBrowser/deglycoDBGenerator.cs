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

        public deglycoDBGenerator(List<string> files, string databasePath)
        {
            this.files = files;
            this.database = databasePath;
        }

        public void Start()
        {
            Dictionary<string, string> db = parseDB(database);

            foreach (string file in files)
            {
                File parsedFile = new File(file, db, Path.GetFileNameWithoutExtension(file));

                printData(parsedFile, file, db);

                printFastaFile(Path.GetDirectoryName(file) + "\\" + Path.GetFileNameWithoutExtension(file) + "_LocalizedPSMs.fasta", 
                    parsedFile.localizedDeglycoPSMsWithMotif, db);
            }            
        }

        public static void printData(File file, string pathout, Dictionary<string, string> db)
        {           

            StreamWriter writer = new StreamWriter(Path.GetDirectoryName(pathout) + "\\" + Path.GetFileNameWithoutExtension(pathout) + "_parsedByGlycoTools.csv");

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

            string uniprotID = "";

            bool isDecoy = true;

            string sequence = "";

            foreach (var fasta in reader.ReadNextFasta())
            {
                if (!fasta.IsDecoy)
                {
                    string desc = fasta.Description.Replace(",", ";");
                    returnDict.Add(">" + desc, fasta.Sequence);
                }                    
            }

            /**
            while (reader.Peek() >= 0)
            {

                string nextLine = reader.ReadLine();

                if (nextLine.StartsWith(">"))
                {

                    if (nextLine.Contains("DECOY"))
                    {
                        isDecoy = true;
                    }
                    else
                    {
                        try
                        {
                            if (!sequence.Equals(""))
                            {
                                returnDict.Add(uniprotID, sequence);
                            }
                            Regex rgx = new Regex(",");
                            uniprotID = rgx.Replace(nextLine, ";");
                            sequence = "";
                            isDecoy = false;
                        }
                        catch (Exception e)
                        {
                            badsequences.Add(sequence);
                            uniprotID = nextLine.Split('|')[1];
                            sequence = "";
                            isDecoy = false;
                        }

                    }
                }
                else
                {
                    if (!isDecoy)
                    {
                        sequence += nextLine;
                    }
                }
            }
    **/

            return returnDict;
        }
    }
}
