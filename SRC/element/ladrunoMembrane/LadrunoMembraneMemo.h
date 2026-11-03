// LADRUNO-HEADER-START
// ==========================================================================
//
//   ▄█          ▄████████ ████████▄     ▄████████ ███    █▄  ███▄▄▄▄    ▄██████▄
//  ███         ███    ███ ███   ▀███   ███    ███ ███    ███ ███▀▀▀██▄ ███    ███
//  ███         ███    ███ ███    ███   ███    ███ ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███  ▄███▄▄▄▄██▀ ███    ███ ███   ███ ███    ███
//  ███       ▀███████████ ███    ███ ▀▀███▀▀▀▀▀   ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███ ▀███████████ ███    ███ ███   ███ ███    ███
//  ███▌    ▄   ███    ███ ███   ▄███   ███    ███ ███    ███ ███   ███ ███    ███
//  █████▄▄██   ███    █▀  ████████▀    ███    ███ ████████▀   ▀█   █▀   ▀██████▀
//  ▀                                   ███    ███
//
//  Ladruno — a research fork of OpenSees
//  Created by:  Nicolas Mora Bowen  ·  Patricio Palacios  ·  José Abell  ·  Guppi
//
// Header auto-stamped by Ladruno_scripts/stamp_headers.py (art: banner_ASCII.txt).
// Do not hand-edit between the markers; edit the script/art and re-run instead.
// ==========================================================================
// LADRUNO-HEADER-END

#ifndef LadrunoMembraneMemo_h
#define LadrunoMembraneMemo_h

// ============================================================================
// LadrunoMembraneMemo.h — WP-222: memo of the pure kernel stages of
// LadrunoMembrane::resolveAndBuildTransfer (the validation mesh, the gravity
// distribution and the diaphragm mesh + LST condensation).
//
// A tower repeats the same panel on every storey: in its own 2D frame the
// inputs of the kernel calls are bit-identical, so the result is too. The key
// is the exact byte image of every input of the stage (doubles by bit pattern,
// so -0.0 and 0.0 are different keys, sizes delimit the lists), and the map
// stores the whole key, so a hit is an exact input match, never a hash
// collision. A stage that throws stores nothing, so it throws again.
//
// The cached values are pure functions of their key, so they can never be
// stale. The memo is still emptied when the last LadrunoMembrane is destroyed
// (wipe(), a rebuilt model), and when it passes a memory budget, so it never
// holds a previous model's panels.
// ============================================================================

#include "LadrunoMembraneKernel.h"

#include <cstddef>
#include <cstring>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace LadrunoMembraneMemo {

class Key {
public:
    void raw(const void *p, std::size_t n) { bytes_.append(static_cast<const char *>(p), n); }
    void d(double x) { raw(&x, sizeof(x)); }
    void i(long long x) { raw(&x, sizeof(x)); }
    void s(const std::string &x) { i(static_cast<long long>(x.size())); bytes_.append(x); }
    void pts(const std::vector<LadrunoMembraneKernel::Point2D> &v) {
        i(static_cast<long long>(v.size()));
        for (const auto &p : v) { d(p.x); d(p.y); }
    }
    void loops(const std::vector<std::vector<LadrunoMembraneKernel::Point2D>> &v) {
        i(static_cast<long long>(v.size()));
        for (const auto &l : v) pts(l);
    }
    void segs(const std::vector<LadrunoMembraneKernel::Segment2D> &v) {
        i(static_cast<long long>(v.size()));
        for (const auto &g : v) { d(g.p1.x); d(g.p1.y); d(g.p2.x); d(g.p2.y); }
    }
    const std::string &bytes() const { return bytes_; }

private:
    std::string bytes_;
};

// one memo per stage; Value must be copy-free to share (held by shared_ptr)
template <class Value>
class Store {
public:
    std::shared_ptr<const Value> find(const Key &k) {
        std::lock_guard<std::mutex> lock(mtx_);
        auto it = map_.find(k.bytes());
        return (it == map_.end()) ? nullptr : it->second;
    }
    void insert(const Key &k, std::shared_ptr<const Value> v, std::size_t valueBytes) {
        std::lock_guard<std::mutex> lock(mtx_);
        std::size_t add = valueBytes + 2 * k.bytes().size() + 64;
        if (bytes_ + add > budget_) {
            map_.clear();
            bytes_ = 0;
        }
        if (map_.emplace(k.bytes(), std::move(v)).second) bytes_ += add;
    }
    void clear() {
        std::lock_guard<std::mutex> lock(mtx_);
        map_.clear();
        bytes_ = 0;
    }
    std::size_t size() {
        std::lock_guard<std::mutex> lock(mtx_);
        return map_.size();
    }

private:
    std::mutex mtx_;
    std::unordered_map<std::string, std::shared_ptr<const Value>> map_;
    std::size_t bytes_ = 0;
    std::size_t budget_ = std::size_t(256) << 20; // per stage
};

template <class T>
inline std::size_t vecBytes(const std::vector<T> &v) { return v.size() * sizeof(T); }

} // namespace LadrunoMembraneMemo

#endif // LadrunoMembraneMemo_h
