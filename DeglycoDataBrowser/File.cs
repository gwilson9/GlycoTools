using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using LumenWorks.Framework.IO.Csv;
using System.IO;
using System.Text.RegularExpressions;

namespace DeglycoDataBrowser
{
    class File
    {
        public string name;
        public List<psms> allPsms;
        public List<psms> psmsWithSequon;
        public List<psms> deglycoPSMs;
        public List<psms> nonLocalizedDeglycoPsms;
        public List<psms> localizedDeglycoPSMsWithMotif;
        public List<psms> TdeglycoPsms;
        public List<psms> SdeglycoPsms;
        public HashSet<psms> uniqueLocalizedPsms;
        public HashSet<string> uniqueLocalizedDeamidatedProteins;


        public File(string pathIn, Dictionary<string, string> db, string name)
        {
            this.name = name;

            this.allPsms = getPsms(pathIn);

            this.deglycoPSMs = parseMods(allPsms);

            this.localizedDeglycoPSMsWithMotif = parseMotifMods(deglycoPSMs, db);

            this.TdeglycoPsms = parseTMods(deglycoPSMs, db);

            this.SdeglycoPsms = parseSMods(deglycoPSMs, db);

            this.psmsWithSequon = getPsmsWithSequon(allPsms);

            this.nonLocalizedDeglycoPsms = getNonLocalizedDeglycoPsms(allPsms);

            this.uniqueLocalizedPsms = new HashSet<psms>(new PsmSequenceComparer());

            this.uniqueLocalizedDeamidatedProteins = new HashSet<string>();

            foreach (psms psm in localizedDeglycoPSMsWithMotif)
            {
                uniqueLocalizedPsms.Add(psm);

                uniqueLocalizedDeamidatedProteins.Add(psm.protein);
            }
        }

        public List<psms> getNonLocalizedDeglycoPsms(List<psms> psms)
        {
            List<psms> returnList = new List<psms>();

            foreach (psms psm in psms)
            {

                string sequence = psm.sequence;

                if (!sequence.Contains("n"))
                {
                    continue;
                }

                bool modSitesWrong = true;

                foreach (Match m in Regex.Matches(sequence, "n"))
                {
                    try
                    {
                        if (sequence[m.Index + 2].Equals('S') || sequence[m.Index + 2].Equals('T'))
                        {
                            modSitesWrong = false;
                            continue;
                        }
                    }
                    catch (Exception e)
                    {

                    }
                }
                if (modSitesWrong)
                {
                    /**
                    foreach(Match m in Regex.Matches(sequence.ToLower(), "n"))
                    {
                        try
                        {
                            if (!sequence[m.Index + 1].Equals('p') && (sequence[m.Index + 2].Equals('s') || sequence[m.Index + 2].Equals('t')))
                            {
                                returnList.Add(psm);
                                continue;
                            }
                        }
                        catch (Exception e)
                        {

                        }
                    }
                    **/
                    foreach (Match m in Regex.Matches(sequence, "N.T"))
                    {
                        returnList.Add(psm);
                        //Console.ReadKey();
                    }

                }
            }

            return returnList;
        }

        public List<psms> getPsmsWithSequon(List<psms> psms)
        {
            List<psms> returnList = new List<psms>();

            foreach (psms psm in psms)
            {
                string sequence = psm.sequence.ToLower();
                foreach (Match m in Regex.Matches(sequence, "n"))
                {
                    try
                    {
                        if (!sequence[m.Index + 1].Equals('p') && (sequence[m.Index + 2].Equals('s') || sequence[m.Index + 2].Equals('t')))
                        {
                            returnList.Add(psm);
                            continue;
                        }
                    }
                    catch (Exception e)
                    {

                    }
                }
            }

            return returnList;
        }

        public List<psms> getPsms(string path)
        {
            List<psms> returnList = new List<psms>();

            using (CsvReader reader = new CsvReader(new StreamReader(path), true))
            {
                while (reader.ReadNextRecord())
                {
                    psms psm = new psms(reader["Spectrum number"], reader["Mass"], reader["Defline"], reader["Peptide"], reader["Charge"], reader["Mods"], reader["Start"], reader["Filename/id"]);
                    returnList.Add(psm);
                }
            }
            return returnList;
        }

        public List<psms> parseMods(List<psms> psms)
        {
            List<psms> returnList = new List<psms>();

            foreach (psms psm in psms)
            {
                if (psm.deglycoMods.Count > 0)
                {
                    returnList.Add(psm);
                }
            }

            return returnList;
        }

        public List<psms> parseMotifMods(List<psms> psms, Dictionary<string, string> db)
        {
            List<psms> returnList = new List<psms>();

            foreach (psms psm in psms)
            {
                foreach (int mod in psm.deglycoMods)
                {
                    string prot = psm.protein.Replace(",", ";");

                    if (mod < db[">" + prot].Length - 1)
                    {
                        if (db[">" + prot][mod + 1].Equals('S') || db[">" + prot][mod + 1].Equals('T'))
                        {
                            returnList.Add(psm);
                        }
                    }
                }
            }
            return returnList;
        }

        public List<psms> parseTMods(List<psms> psms, Dictionary<string, string> db)
        {
            List<psms> returnList = new List<psms>();

            foreach (psms psm in psms)
            {
                foreach (int mod in psm.deglycoMods)
                {
                    string prot = psm.protein.Replace(",", ";");

                    if (mod < db[">" + prot].Length - 1)
                    {
                        if (db[">" + prot][mod + 1].Equals('T'))
                        {
                            returnList.Add(psm);
                        }
                    }

                }
            }
            return returnList;
        }

        public List<psms> parseSMods(List<psms> psms, Dictionary<string, string> db)
        {
            List<psms> returnList = new List<psms>();

            foreach (psms psm in psms)
            {
                foreach (int mod in psm.deglycoMods)
                {
                    string prot = psm.protein.Replace(",", ";");

                    if (mod < db[">" + prot].Length - 1)
                    {
                        if (db[">" + prot][mod + 1].Equals('S'))
                        {
                            returnList.Add(psm);
                        }
                    }
                }
            }
            return returnList;
        }
    }
}
