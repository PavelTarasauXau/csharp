using System;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations;
using System.IO;
using System.Linq.Expressions;

namespace GeneticsProject
{
    public struct GeneticData
    {
        public string name; //protein name
        public string organism;
        public string formula; //formula
    }

    class Program
    {
        static List<GeneticData> data = new List<GeneticData>();
        static int count = 1;
        static string GetFormula(string proteinName)
        {
            foreach (GeneticData item in data)
            {
                if (item.name.Equals(proteinName)) return item.formula;
            }
            return null;
        }
        static void ReadGeneticData(string filename)
        {
            StreamReader reader = new StreamReader(filename);
            while (!reader.EndOfStream)
            {
                string line = reader.ReadLine();
                string[] fragments = line.Split('\t');
                GeneticData protein;
                protein.name = fragments[0];
                protein.organism = fragments[1];
                protein.formula = fragments[2];
                data.Add(protein);
                count++;
            }
            reader.Close();
        }
        static void ReadHandleCommands(string filename)
        {
            StreamWriter genedataWr = new StreamWriter("genedata.txt"); 
            StreamReader reader = new StreamReader(filename);
            genedataWr.WriteLine("Tarasov Pavel");
            genedataWr.WriteLine("Genetic Searching");
            genedataWr.WriteLine("================================================");
            int counter = 0;
            while (!reader.EndOfStream)
            {
                string line = reader.ReadLine(); counter++;
                string[] command = line.Split('\t');
                if (command[0].Equals("search"))
                {
                    //001   search  SIIK
                    genedataWr.WriteLine($"{counter.ToString("D3")}   {"search"}   {Decoding(command[1])}");
                    int index = Search(command[1]);
                    if (index != -1)
                        genedataWr.WriteLine($"organism: {data[index].organism}    protein: {data[index].name}");
                    else
                        genedataWr.WriteLine("NOT FOUND");
                    genedataWr.WriteLine("================================================");
                }
                if (command[0].Equals("diff"))
                {
                    genedataWr.WriteLine($"{counter.ToString("D3")}   {"diff"}   {command[1]}    {command[2]}");

                    string[] proteins = new string[2];
                      
                    int ind = 0;
                    foreach (GeneticData elem in data)
                    {
                        if (elem.name == command[1]|| elem.name == command[2])
                        {
                            proteins[ind]=elem.formula;
                            ind++;
                            

                        }
                        
                    }
                    if ((ind != 2 && command[1] != command[2] )|| (command[1] == command[2]&&ind!=1))
                    {
                        genedataWr.WriteLine("amino-acids difference: MISSING");                   
                    }
                    else
                    {
                        int result = 0;
                        if (command[1] != command[2]) result = Diff(Decoding(proteins[0]), Decoding(proteins[1]));
                        

                        genedataWr.WriteLine("amino-acids difference: "+result);

                    }
                    genedataWr.WriteLine("================================================");



                }
                if (command[0].Equals("mode"))
                {
                    genedataWr.WriteLine($"{counter.ToString("D3")}   {"mode"}   {Decoding(command[1])}");
                    genedataWr.WriteLine("amino-acid occurs:  ");
                    string s = "";
                    int ind = 0;
                    foreach (GeneticData elem in data)
                    {
                        if (elem.name == command[1])
                        {
                            s = elem.formula;
                            ind++;
                            break;


                        }

                    }
                    if (ind != 1)
                    {
                        genedataWr.WriteLine("amino-acids difference: MISSING");
                    }
                    else
                    {
                        string result = Mode(Decoding(s));
                        genedataWr.WriteLine(result);

                    }

                    genedataWr.WriteLine("================================================");
                }
            }
            reader.Close();
            genedataWr.Close();
        }
        static bool IsValid(string formula)
        {
            List<char> letters = new List<char>() { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' };
            foreach (char ch in formula)
            {
                if (!letters.Contains(ch)) return false;
            }
            return true;
        }
        static string Encoding(string formula)
        {
            string encoded = String.Empty;
            for (int i = 0; i < formula.Length; i++)
            {
                char ch = formula[i];
                int count = 1;
                while (i < formula.Length - 1 && formula[i + 1] == ch)
                {
                    count++;
                    i++;
                }
                if (count > 2) encoded = encoded + count + ch;
                if (count == 1) encoded = encoded + ch;
                if (count == 2) encoded = encoded + ch + ch;

            }
            return encoded;


        }
        static string Decoding(string formula)
        {
            string decoded = String.Empty;
            for (int i = 0; i < formula.Length; i++)
            {
                if (char.IsDigit(formula[i]))
                {
                    char letter = formula[i + 1];
                    int conversion = formula[i] - '0';
                    for (int j = 0; j < conversion - 1; j++) decoded = decoded + letter;
                }
                else decoded = decoded + formula[i];
            }
            return decoded;
        }
        static int Search(string amino_acid)
        {
            //       FKIII                FK3I
            string decoded = Decoding(amino_acid);
            for (int i = 0; i < data.Count; i++)
            {
                if (data[i].formula.Contains(decoded)) return i;
            }
            return -1;
        }
        static int Diff(string protein1, string protein2)
        {
            int c = 0; 
            int minLength = Math.Min(protein1.Length, protein2.Length);
            for (int i = 0; i < minLength; i++)
            {
                if (protein1[i] != protein2[i])
                {
                    c++;
                    
                }
                
            }
            return c+Math.Abs(protein1.Length-protein2.Length);
        }
        static string Mode(string protein)
        {
            List<char> ltrs = new List<char>() { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' };
            Dictionary<char, int> freq = new Dictionary<char, int>();
            string s = "";
            foreach (char ltr in ltrs)
            {
                freq[ltr] = 0;
            }
            foreach (char c in protein)
            {
                if (freq.ContainsKey(c))
                {
                    freq[c]++;
                }
            }
            char mostFreqChar = '\0';
            int maxFreq = 0;

            foreach (char ltr in ltrs)
            {
                if (freq[ltr] > maxFreq)
                {
                    mostFreqChar = ltr;
                    maxFreq = freq[ltr];
                }
            }
            s = $"{mostFreqChar}\t{maxFreq}";
            return s;

        }
        static void Main(string[] args)
        {

            
            ReadGeneticData("sequences\\sequences.2.txt");
            
            ReadHandleCommands("commands\\commands.2.txt");
            
            
        }
    }
}