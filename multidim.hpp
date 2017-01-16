#ifndef CLASSIC_HPP
#define CLASSIC_HPP
#include <iostream>
#include <array>
#include <vector>
#include <utility>
#include <unordered_map>
#include "common.h"
#include "protein.hpp"

using namespace std;

template <size_t D>
struct MultiDim {
    template<typename T>
    class Vec : public std::array<T, D> {
    public:
        constexpr double operator- (const Vec &v) const {
            double sum = 0;
            for (size_t i=0; i<D; i++) {
                double d = (*this)[i] - v[i];
                sum += d*d;
            }
            return sqrt(sum);
        }

        constexpr bool operator== (const Vec &v) const {
            bool ok = true;
            for (size_t i=0; i<D; i++)
                ok &= ((*this)[i] == v[i]);
            return ok;
        }

        constexpr bool operator!= (const Vec &v) const {
            return !((*this) == v);
        }

        explicit constexpr operator double() const {
            double sum = 0;
            for (size_t i=0; i<D; i++) {
                double d = (*this)[i];
                sum += d*d;
            }
            return sqrt(sum);
        }

        explicit constexpr operator int() const {
            int h = 0;
            for (size_t i=0; i<D; i++)
                h = h ^ (*this)[i];
            return h;
        }

        template <typename ...E>
        Vec(E&& ...l) : std::array<T,D>{{std::forward(l)...}} {}

        Vec (const nlohmann::json &j) : std::array<T,D>(j) {}

//        constexpr operator nlohmann::json() const {
//            return "";
//        }
    };

public:
    // we use 3 coordinates proteins (id, enh, inh)
    using ProteinInternalCoord_t = int;
    using ProteinCoord_t = Vec<ProteinInternalCoord_t>;
    static constexpr ProteinInternalCoord_t PROTEIN_LOWER = 0;
    static constexpr ProteinInternalCoord_t PROTEIN_UPPER = 32;
    static constexpr double PROTEIN_MAX_AMPLITUDE = (PROTEIN_UPPER + PROTEIN_LOWER) / 2.;
    static constexpr double PROTEIN_MAX_DISTANCE = double(ProteinCoord_t((PROTEIN_LOWER + PROTEIN_LOWER) / 2.));
    using Protein_t = Protein<3, ProteinCoord_t, PROTEIN_LOWER, PROTEIN_UPPER>;

    // we need 2 parameters (beta, alpha)
    static constexpr unsigned int nbParams = 2;
    // and we produce 2 dimensional signatures (enhnance, inhibit)
    static constexpr unsigned int nbSignatureParams = 2;

    static const array<pair<double, double>, nbParams> paramsLimits() {
        return {{{0.5, 2.0}, {0.5, 2.0}}};
    }

    // helpers for proteins coords
    static inline ProteinCoord_t& getId(Protein_t& p) { return p.coords[0]; }
    static inline ProteinCoord_t& getEnh(Protein_t& p) { return p.coords[1]; }
    static inline ProteinCoord_t& getInh(Protein_t& p) { return p.coords[2]; }

    // aliases for ProteinType
    static constexpr ProteinType pinput = ProteinType::input;
    static constexpr ProteinType pregul = ProteinType::regul;
    static constexpr ProteinType poutput = ProteinType::output;

    double maxEnhance = 0.0, maxInhibit = 0.0;

    MultiDim() {}

    template <typename GRN> void updateSignatures(GRN& grn) {
        grn.signatures.clear();
        grn.signatures.resize(grn.actualProteins.size());
        for (size_t i = 0; i < grn.actualProteins.size(); ++i) {
            grn.signatures[i].resize(grn.actualProteins.size());
            for (size_t j = 0; j < grn.actualProteins.size(); ++j) {
                auto& p0 = grn.actualProteins[i];
                auto& p1 = grn.actualProteins[j];
                grn.signatures[i][j] = {
                    {similarity(getEnh(p0), getId(p1)),
                     similarity(getInh(p0), getId(p1))}
                };
                if (grn.signatures[i][j][0] > maxEnhance) maxEnhance = grn.signatures[i][j][0];
                if (grn.signatures[i][j][1] > maxInhibit) maxInhibit = grn.signatures[i][j][1];
            }
        }
        // std::cerr << "maxEnh = " << maxEnhance << ", maxInh = " << maxInhibit << std::endl;
        for (size_t i = 0; i < grn.actualProteins.size(); ++i) {
            for (size_t j = 0; j < grn.actualProteins.size(); ++j) {
                grn.signatures[i][j] = {
                    {exp(grn.params[0] * grn.signatures[i][j][0] - maxEnhance),
                     exp(grn.params[0] * grn.signatures[i][j][1] - maxInhibit)}};
            }
        }
    }

    static constexpr double similarity (const ProteinCoord_t &c0, const ProteinCoord_t &c1) {
        return static_cast<double>(
            PROTEIN_MAX_AMPLITUDE - std::min(fabs(c0 - c1), fabs(c1 - c0))
        );
    }

    template <typename GRN> void step(GRN& grn, unsigned int nbSteps) {
        for (auto s = 0u; s < nbSteps; ++s) {
            std::vector<double> nextProteins;  // only reguls & outputs concentrations
            nextProteins.reserve(grn.getNbProteins() - grn.getProteinSize(ProteinType::input));
            const auto firstOutputId = grn.getFirstOutputIndex();
            for (size_t j = grn.getFirstRegulIndex(); j < grn.getNbProteins(); ++j) {
                double enh = 0.0, inh = 0.0;
                for (size_t k = 0; k < firstOutputId; ++k) {
                    enh += grn.actualProteins[k].c * grn.signatures[k][j][0];
                    inh += grn.actualProteins[k].c * grn.signatures[k][j][1];
                }
                nextProteins.push_back(
                            max(0.0, grn.actualProteins[j].c +
                                (grn.params[1] / static_cast<double>(grn.getNbProteins())) *
                            (enh - inh)));
            }
            // Normalizing regul & output proteins concentrations
            double sumConcentration = 0.0;
            for (auto i : nextProteins) {
                sumConcentration += i;
            }
            if (sumConcentration > 0) {
                for (auto& i : nextProteins) {
                    i /= sumConcentration;
                }
            }
            auto firstRegulIndex = grn.getFirstRegulIndex();
            for (size_t i = firstRegulIndex; i < grn.getNbProteins(); ++i) {
                grn.actualProteins[i].c = nextProteins[i - firstRegulIndex];
            }
        }
    }

    static constexpr ProteinCoord_t getRandomCoord (void) {
        using dist_t = std::conditional<
            std::is_integral<ProteinInternalCoord_t>::value,
            std::uniform_int_distribution<ProteinInternalCoord_t>,
            std::uniform_real_distribution<ProteinInternalCoord_t>
        >::type;

        ProteinCoord_t vec;
        dist_t dist (PROTEIN_LOWER, PROTEIN_UPPER);
        for (auto &val: vec)
            val = dist(grnRand);
        return vec;
    }

};
#endif
