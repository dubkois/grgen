#ifndef _MULTI_DIMENSIONS_HPP_
#define _MULTI_DIMENSIONS_HPP_

#include <cstring>

#include "common.h"
#include "protein.hpp"

template <typename T, size_t D, int LOWER=-1, int UPPER=1>
struct MultiDim {
    using json = nlohmann::json;

    class Vec {
        static constexpr T AMPLITUDE = UPPER - LOWER;
//        static constexpr T HALF_AMPLITUDE = AMPLITUDE / 2.;

        static constexpr double CORNERS_DISTANCE = sqrt(D * AMPLITUDE*AMPLITUDE);

        std::array<T,D> _data;

    public:
        static constexpr double MAX_DISTANCE = CORNERS_DISTANCE / 2.;

        constexpr double operator- (const Vec &v) const {
            double sum = 0;
            for (size_t i=0; i<D; i++) {
                double d = _data[i] - v._data[i];
                sum += d*d;
            }
            sum = sqrt(sum);
            return std::min(sum, CORNERS_DISTANCE - sum);
        }

        explicit constexpr operator double() const {
            double sum = 0;
            for (size_t i=0; i<D; i++) {
                double d = _data[i];
                sum += d*d;
            }
            return sqrt(sum);
        }

        explicit constexpr operator int() const {
            int h = 0;
            for (size_t i=0; i<D; i++)
                h = h ^ _data[i];
            return h;
        }

        constexpr Vec (T val=T((LOWER+UPPER)/2.))
            : _data ({val}) {}

        Vec (const json &j) {
            auto vcoords = j.get<std::vector<T>>();
            std::copy(vcoords.begin(), vcoords.end(), _data.begin());
        }

        explicit constexpr operator json() const {
            return json(_data);
        }

        friend std::ostream& operator<< (std::ostream &os, const Vec &v) {
            std::ostream_iterator<T> output(os, " ");
            os << "(";
            std::copy(v._data.begin(), v._data.end(), output);
            return os << ")";
        }

        static constexpr Vec getRandomCoord (void) {
            using dist_t = typename std::conditional<
                std::is_integral<T>::value,
                std::uniform_int_distribution<T>,
                std::uniform_real_distribution<T>
            >::type;

            Vec vec;
            dist_t dist (LOWER, UPPER);
            for (auto &val: vec._data)
                val = dist(grnRand);
            return vec;
        }
    };

public:
    // we use 3 coordinates proteins (id, enh, inh)
    using ProteinInternalCoord_t = T;
    using ProteinCoord_t = Vec;
    using Protein_t = Protein<3, ProteinCoord_t, LOWER, UPPER>;

    // we need 2 parameters (beta, alpha)
    static constexpr unsigned int nbParams = 2;
    // and we produce 2 dimensional signatures (enhnance, inhibit)
    static constexpr unsigned int nbSignatureParams = 2;

    static const std::array<std::pair<double, double>, nbParams> paramsLimits() {
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
//        std::cout << "sim(" << c0 << ", " << c1 << ") = "
//                  << Vec::MAX_DISTANCE << " - ("
//                  << c0 - c1
//                  << ") = "
//                  << Vec::MAX_DISTANCE - (c0 - c1)
//                  << std::endl;
        return static_cast<double>(
            Vec::MAX_DISTANCE - (c0 - c1)
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
};

#endif
