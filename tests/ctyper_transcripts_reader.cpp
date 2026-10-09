// Optional integration harness. Compile against an unchanged Ctyper checkout;
// see docs/transcripts.md. This exercises its real index, tree and row readers.
#include "PriorData.hpp"
#include <cassert>

std::string matrix_version = MATRIX_VERSION_V2;
MatrixLayout matrix_layout = MatrixLayout::V2_FIXED;
bool optioncorr = false;

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    PriorData reader(argv[1], 1);
    const size_t chunks = reader.LoadIndex();
    for (size_t i = 0; i < chunks; ++i) {
        PriorChunk* chunk = reader.getChunkData(i);
        assert(chunk->genenum > 0);
        assert(chunk->pathnames.size() == 1);  // native no-path sentinel
        assert(chunk->numgroups == chunk->genenum);
        size_t leaves = 0;
        for (size_t n = 0; n < 2 * chunk->genenum - 1; ++n) {
            const node& item = chunk->phylo_tree[n];
            assert(item.numchildren == 0 || item.numchildren == 2);
            if (!item.numchildren) assert(item.index == leaves++);
        }
        assert(leaves == chunk->genenum);
        std::cout << chunk->prefix << '\t' << chunk->genenum << '\t' << chunk->kmervec_size << '\n';
        uint16* row = chunk->kmer_matrix;
        std::vector<unsigned long long> totals(chunk->genenum, 0);
        for (size_t k = 0; k < chunk->kmervec_size; ++k) {
            assert(row[0] == '+');
            assert(row[2] == 0 && row[3] == 0 && row[4] == 0);
            assert(row[5] == (100 << 8) + (100 << 1));
            std::cout << "ROW";
            for (size_t j = 0; j < row[1]; ++j) {
                const auto allele = row[FIXCOL + j];
                assert(allele < chunk->genenum);
                totals[allele]++;
                std::cout << '\t' << allele;
            }
            std::cout << '\n';
            row += FIXCOL + row[1];
        }
        for (size_t j = 0; j < totals.size(); ++j) assert(totals[j] == chunk->gene_kmercounts[j]);
        reader.FinishChunk(chunk);
    }
    return 0;
}
