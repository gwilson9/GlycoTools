using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DeglycoDataBrowser
{
    class HighlightEventArgs : EventArgs
    {
        public string rawFile { get; private set; }

        public string resultsFile { get; private set; }

        public HighlightEventArgs(string resultsFile)
        {
            this.resultsFile = resultsFile;
        }
    }
}
