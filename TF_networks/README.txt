All network information can be extracted by running: 'tar -xzf networks.tgz'

Each directory (cell type, n=41) includes a file (genes-regulate-genes.txt)
containing all of the TF-to-TF regulatory interactions mapped within that cell
type.

Regulatory interactions were identified in each cell type by scanning the
proximal DHSs (+/-5kb from the canonical transcriptional start site) of each
transcription factor gene DNaseI footprints corresponding to the recognition
sequences of known TFs. On average, 11,193 regulatory interactions were
identified in each cell type.

The file genes-regulate-genes.txt contains two columns of gene names. The
gene in the first column is bound (ie contains a DNaseI footprint in proximal
regulatory DNA) corresponding to the second column's TF gene product.


Note: regulatory direction: Column 2 --> Column 1
