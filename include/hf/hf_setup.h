#include "hf/chem.h"
#include <cstdint>
#include <optional>

namespace hf
{
    auto cgf_overlap(const contracted_gaussian_functions& cgf1, const contracted_gaussian_functions& cgf2) -> double;
    auto cgf_kinetic(const contracted_gaussian_functions& cgf1, const contracted_gaussian_functions& cgf2) -> double;
    auto cgf_nuclear(contracted_gaussian_functions& cgf1, contracted_gaussian_functions& cgf2, const atom& a) -> double;
    auto cgf_repulsion(const contracted_gaussian_functions& cgf1, const contracted_gaussian_functions& cgf2,
        const contracted_gaussian_functions& cgf3, const contracted_gaussian_functions& cgf4) -> double;
    constexpr auto te_index(uint32_t i, uint32_t j, uint32_t k, uint32_t l) -> uint32_t
    {
        if (i < j)
        {
            std::swap(i, j);
        }
        if (k < l)
        {
            std::swap(k, l);
        }
        uint32_t ij = i * (i + 1) / 2 + j;
        uint32_t kl = k * (k + 1) / 2 + l;
        if (ij < kl)
        {
            std::swap(ij, kl);
        }
        return ij * (ij + 1) / 2 + kl;
    }
} // namespace hf
