#pragma once

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <ostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Grid diagrams / arc presentations.
//
// A GRID DIAGRAM of size n is an n x n array of cells carrying exactly one X
// and one O in every row and every column.  Join the two markers of each row
// by a horizontal segment and the two markers of each column by a vertical
// segment; declare every vertical to pass OVER every horizontal.  The result
// is a knot diagram.
//
// The same object read the other way is an ARC PRESENTATION: the columns are
// the n arcs (one per half-plane, or "page", around a binding axis), and the
// rows are the n binding points where consecutive arcs meet on the axis.  The
// horizontal segments are the pieces of the axis and therefore run behind
// everything -- that is exactly the "verticals over horizontals" rule.
//
// The minimal n over all grid diagrams of a knot K is the ARC INDEX alpha(K)
// (equivalently the grid number).  Nothing here is geometric: a grid diagram
// is a pair of permutations, which is the same kind of data as a petal
// permutation.  Indeed:
//
//   PETAL PERMUTATION  ==>  ARC PRESENTATION  ==>  GRID DIAGRAM
//
// A petal projection with n petals has the knot passing through the central
// ubercrossing n times, at heights pi_1, ..., pi_n in traversal order, with
// the t-th petal joining the passage at height pi_t to the one at height
// pi_{t+1}.  Read the central axis as the binding axis: the n passages are
// the n binding points (row pi_t) and the n petals are the n arcs (one
// column each).  Petals are laid out around the centre with the traversal
// step s = (n+1)/2, so petal t sits in angular slot (t*s) mod n -- that is
// the column index.  Hence
//
//      alpha(K) <= n   for every knot K from a length-n petal permutation,
//
// and in particular alpha(K) <= p(K), the petal number.
//
// Representation below: for each column j, xrow_[j] and orow_[j] are the rows
// of its X and O.  Both are permutations of 0..n-1.  Orientation convention:
// verticals run X -> O, horizontals run O -> X, giving a coherently oriented
// closed curve.
// ---------------------------------------------------------------------------
class GridDiagram {
public:
    GridDiagram(std::vector<int> xrow, std::vector<int> orow)
        : xrow_(std::move(xrow)), orow_(std::move(orow)) {
        validate();
    }

    int size() const { return static_cast<int>(xrow_.size()); }

    int xRowOfCol(int j) const { return xrow_[j]; }
    int oRowOfCol(int j) const { return orow_[j]; }

    // Column carrying the X (resp. O) of row i -- the inverse permutations.
    int xColOfRow(int i) const { return xcol_()[i]; }
    int oColOfRow(int i) const { return ocol_()[i]; }

    // -----------------------------------------------------------------
    // Petal permutation -> grid diagram.
    //
    // perm is 1-based, of odd length n; perm[t] is the height of the t-th
    // passage through the centre in traversal order.  Arc t joins heights
    // perm[t] and perm[t+1] and lives in column (t*s) mod n, s = (n+1)/2.
    // Orienting the knot along the traversal puts the X at row perm[t]-1
    // and the O at row perm[t+1]-1.
    // -----------------------------------------------------------------
    static GridDiagram fromPetalPermutation(const std::vector<int>& perm) {
        const int n = static_cast<int>(perm.size());
        if (n < 3 || n % 2 == 0)
            throw std::invalid_argument("petal permutation length must be odd and >= 3");
        const int s = (n + 1) / 2;
        std::vector<int> xr(n), orr(n);
        for (int t = 0; t < n; ++t) {
            const int col = (t * s) % n;
            xr[col]  = perm[t] - 1;
            orr[col] = perm[(t + 1) % n] - 1;
        }
        return GridDiagram(std::move(xr), std::move(orr));
    }

    // Smallest legal grid: the 2x2 unknot.
    static GridDiagram unknot() { return GridDiagram({0, 1}, {1, 0}); }

    // -----------------------------------------------------------------
    // Diagram statistics
    // -----------------------------------------------------------------

    // Crossings of the underlying diagram: one for every (row, column) pair
    // whose horizontal and vertical segments strictly straddle each other.
    int crossingCount() const {
        const int n = size();
        int count = 0;
        for (int i = 0; i < n; ++i) {
            const int ha = std::min(xColOfRow(i), oColOfRow(i));
            const int hb = std::max(xColOfRow(i), oColOfRow(i));
            for (int j = ha + 1; j < hb; ++j) {
                const int va = std::min(xrow_[j], orow_[j]);
                const int vb = std::max(xrow_[j], orow_[j]);
                if (va < i && i < vb) ++count;
            }
        }
        return count;
    }

    // Number of link components; 1 means we are looking at a knot.
    int componentCount() const {
        const int n = size();
        std::vector<char> seen(n, 0);
        int comps = 0;
        for (int start = 0; start < n; ++start) {
            if (seen[start]) continue;
            ++comps;
            int j = start;
            do {
                seen[j] = 1;
                j = xColOfRow(orow_[j]);   // vertical X->O, then horizontal O->X
            } while (j != start);
        }
        return comps;
    }

    // -----------------------------------------------------------------
    // Cromwell moves
    // -----------------------------------------------------------------

    // Cyclic permutation (translation on the torus).
    void rotateColumns(int k) {
        const int n = size();
        k = ((k % n) + n) % n;
        if (k == 0) return;
        std::vector<int> nx(n), no(n);
        for (int j = 0; j < n; ++j) {
            nx[(j + k) % n] = xrow_[j];
            no[(j + k) % n] = orow_[j];
        }
        xrow_ = std::move(nx); orow_ = std::move(no);
        invalidate();
    }

    void rotateRows(int k) {
        const int n = size();
        k = ((k % n) + n) % n;
        if (k == 0) return;
        for (int j = 0; j < n; ++j) {
            xrow_[j] = (xrow_[j] + k) % n;
            orow_[j] = (orow_[j] + k) % n;
        }
        invalidate();
    }

    // Commutation of columns j and j+1: legal when the two vertical segments'
    // row intervals are disjoint or nested (not interleaved).
    bool commuteColumns(int j) {
        const int n = size();
        const int k = (j + 1) % n;
        if (!intervalsCommute(xrow_[j], orow_[j], xrow_[k], orow_[k])) return false;
        std::swap(xrow_[j], xrow_[k]);
        std::swap(orow_[j], orow_[k]);
        invalidate();
        return true;
    }

    // Commutation of rows i and i+1, by the same test on column intervals.
    bool commuteRows(int i) {
        const int n = size();
        const int k = (i + 1) % n;
        if (!intervalsCommute(xColOfRow(i), oColOfRow(i),
                              xColOfRow(k), oColOfRow(k))) return false;
        for (int j = 0; j < n; ++j) {
            if (xrow_[j] == i)      xrow_[j] = k;
            else if (xrow_[j] == k) xrow_[j] = i;
            if (orow_[j] == i)      orow_[j] = k;
            else if (orow_[j] == k) orow_[j] = i;
        }
        invalidate();
        return true;
    }

    // Destabilization: find a 2x2 block (cyclically adjacent rows and
    // columns) holding exactly three markers.  The marker diagonally
    // opposite the empty cell is the "corner"; deleting its row and column
    // and dropping a marker of the shared type into the empty cell reduces
    // the grid size by one.  Returns false if no such block exists.
    bool tryDestabilize() {
        const int n = size();
        if (n <= 2) return false;
        for (int i = 0; i < n; ++i) {
            const int i2 = (i + 1) % n;
            for (int j = 0; j < n; ++j) {
                const int j2 = (j + 1) % n;
                const int m[4] = { markerAt(i, j),  markerAt(i, j2),
                                   markerAt(i2, j), markerAt(i2, j2) };
                int filled = 0, emptyIdx = -1;
                for (int q = 0; q < 4; ++q) {
                    if (m[q]) ++filled; else emptyIdx = q;
                }
                if (filled != 3) continue;

                const int rows[2] = { i, i2 };
                const int cols[2] = { j, j2 };
                const int er = rows[emptyIdx / 2], ec = cols[emptyIdx % 2];
                const int cr = rows[1 - emptyIdx / 2], cc = cols[1 - emptyIdx % 2];
                // The two non-corner markers share a type; that type moves
                // into the empty cell.
                const int type = markerAt(er, cc);
                doDestabilize(cr, cc, er, ec, type);
                return true;
            }
        }
        return false;
    }

    // Stabilization: the inverse of the move above, and the only Cromwell
    // move that increases the grid size. Pick the marker of type `type` in
    // column j; insert a new row just above its row and a new column just
    // right of j, and expand the marker into a 2x2 block carrying three
    // markers. Needed because monotonic simplification is not complete: a
    // grid can be stuck at a size above the arc index with no
    // destabilization available in any commutation class reachable from it.
    void stabilize(int j, bool onX) {
        const int n = size();
        const int r = onX ? xrow_[j] : orow_[j];
        const int rNew = r + 1;                  // inserted row index
        const int cNew = j + 1;                  // inserted column index

        auto liftRow = [r](int x) { return x > r ? x + 1 : x; };

        std::vector<int> nx(n + 1), no(n + 1);
        for (int c = 0; c < n; ++c) {
            const int dst = (c < cNew) ? c : c + 1;
            nx[dst] = liftRow(xrow_[c]);
            no[dst] = liftRow(orow_[c]);
        }
        // Old column j keeps its partner marker and takes the new row for the
        // marker we stabilized at; the new column carries the other two.
        if (onX) { nx[j] = rNew;  nx[cNew] = r;  no[cNew] = rNew; }
        else     { no[j] = rNew;  no[cNew] = r;  nx[cNew] = rNew; }

        xrow_ = std::move(nx); orow_ = std::move(no);
        invalidate();
        validate();
    }

    // -----------------------------------------------------------------
    // Arc-index search.
    //
    // Destabilize greedily; when stuck, wander with random commutations and
    // try again.  Commutations preserve the grid size and destabilizations
    // shrink it, so the search is monotone in the best size seen.  Returns
    // that best size -- an UPPER BOUND for the arc index, and in practice
    // the exact value for the small knots this pipeline produces.
    // -----------------------------------------------------------------
    int simplify(std::uint64_t seed = 12345u, int steps = 40000, int slack = 2) {
        std::mt19937_64 rng(seed);
        greedyDestabilize();
        GridDiagram best = *this;

        for (int step = 0; step < steps && size() > 2; ++step) {
            const int n = size();
            std::uniform_int_distribution<int> pick(0, n - 1);
            const int roll = static_cast<int>(rng() % 100);

            if (roll < 80) {
                // Commutations: size-preserving, they shuffle the diagram
                // within its commutation class looking for a corner.
                if (rng() & 1u) commuteColumns(pick(rng));
                else            commuteRows(pick(rng));
            } else if (size() < best.size() + slack) {
                // Go uphill occasionally. A stabilization followed by
                // commutations can expose a destabilization that was not
                // available before, which is how the search gets past a
                // plateau.
                stabilize(pick(rng), (rng() & 1u) != 0);
                for (int k = 0, kn = size(); k < kn; ++k) {
                    std::uniform_int_distribution<int> p2(0, size() - 1);
                    if (rng() & 1u) commuteColumns(p2(rng));
                    else            commuteRows(p2(rng));
                }
            }

            greedyDestabilize();
            if (size() < best.size()) { best = *this; step = 0; }
            else if (size() > best.size() + slack) { *this = best; }
        }

        if (best.size() < size()) *this = best;
        return size();
    }

    void greedyDestabilize() { while (tryDestabilize()) {} }

    // -----------------------------------------------------------------
    // Output
    // -----------------------------------------------------------------

    // ASCII picture, row 0 at the bottom.
    std::string toAscii() const {
        const int n = size();
        std::ostringstream out;
        for (int i = n - 1; i >= 0; --i) {
            for (int j = 0; j < n; ++j) {
                const int m = markerAt(i, j);
                out << (m == kX ? 'X' : m == kO ? 'O' : '.');
                if (j + 1 < n) out << ' ';
            }
            out << '\n';
        }
        return out.str();
    }

    // One line per column: "col xrow orow".
    std::string toPermutationPair() const {
        const int n = size();
        std::ostringstream out;
        out << "X = [";
        for (int j = 0; j < n; ++j) out << (j ? " " : "") << xrow_[j] + 1;
        out << "], O = [";
        for (int j = 0; j < n; ++j) out << (j ? " " : "") << orow_[j] + 1;
        out << "]";
        return out.str();
    }

    // Closed polygon in R^3 realising the diagram, in the same TSV format
    // petal_star_clean emits, so it can be piped straight into
    // knoodlesimplify / knoodleidentify.
    //
    // Horizontals sit at z = -1, verticals at z = +1, so verticals really do
    // pass over horizontals.  Each corner is cut by a short diagonal of
    // length eps carrying the height change, which keeps every projected
    // edge non-degenerate.
    void writePolygonTSV(std::ostream& out, double eps = 0.25) const {
        out << std::fixed << std::setprecision(10);
        int j = 0;
        do {
            const int r1 = xrow_[j], r2 = orow_[j];
            const int prevCol = oColOfRow(r1);
            const int nextCol = xColOfRow(r2);
            const double sxIn  = sgn(j - prevCol);
            const double sxOut = sgn(nextCol - j);
            const double sy    = sgn(r2 - r1);

            emit(out, j - eps * sxIn,  r1,              -1.0);  // end of horizontal
            emit(out, j,               r1 + eps * sy,   +1.0);  // start of vertical
            emit(out, j,               r2 - eps * sy,   +1.0);  // end of vertical
            emit(out, j + eps * sxOut, r2,              -1.0);  // start of horizontal
            j = nextCol;
        } while (j != 0);
    }

private:
    static constexpr int kNone = 0, kX = 1, kO = 2;

    std::vector<int> xrow_, orow_;
    mutable std::vector<int> xcol_cache_, ocol_cache_;

    static double sgn(int v) { return v > 0 ? 1.0 : (v < 0 ? -1.0 : 0.0); }

    static void emit(std::ostream& out, double x, double y, double z) {
        out << x << '\t' << y << '\t' << z << '\n';
    }

    void invalidate() const { xcol_cache_.clear(); ocol_cache_.clear(); }

    const std::vector<int>& xcol_() const {
        if (xcol_cache_.empty()) {
            xcol_cache_.assign(size(), 0);
            for (int j = 0; j < size(); ++j) xcol_cache_[xrow_[j]] = j;
        }
        return xcol_cache_;
    }

    const std::vector<int>& ocol_() const {
        if (ocol_cache_.empty()) {
            ocol_cache_.assign(size(), 0);
            for (int j = 0; j < size(); ++j) ocol_cache_[orow_[j]] = j;
        }
        return ocol_cache_;
    }

    int markerAt(int i, int j) const {
        if (xrow_[j] == i) return kX;
        if (orow_[j] == i) return kO;
        return kNone;
    }

    // Non-interleaved test used by both commutation moves.
    static bool intervalsCommute(int a1, int a2, int b1, int b2) {
        const int alo = std::min(a1, a2), ahi = std::max(a1, a2);
        const int blo = std::min(b1, b2), bhi = std::max(b1, b2);
        const bool disjoint = (ahi < blo) || (bhi < alo);
        const bool nested   = (alo < blo && bhi < ahi) || (blo < alo && ahi < bhi);
        return disjoint || nested;
    }

    // Delete row cr and column cc, after placing a marker of the given type
    // in cell (er, ec).
    void doDestabilize(int cr, int cc, int er, int ec, int type) {
        const int n = size();
        std::vector<int> nx, no;
        nx.reserve(n - 1); no.reserve(n - 1);
        auto shiftRow = [cr](int r) { return r > cr ? r - 1 : r; };
        for (int j = 0; j < n; ++j) {
            if (j == cc) continue;
            int xr = xrow_[j], orr = orow_[j];
            if (j == ec) { if (type == kX) xr = er; else orr = er; }
            nx.push_back(shiftRow(xr));
            no.push_back(shiftRow(orr));
        }
        xrow_ = std::move(nx); orow_ = std::move(no);
        invalidate();
        validate();
    }

    void validate() const {
        const int n = size();
        if (n < 2) throw std::invalid_argument("grid size must be at least 2");
        if (static_cast<int>(orow_.size()) != n)
            throw std::invalid_argument("X and O row lists must have equal length");
        std::vector<char> seenX(n, 0), seenO(n, 0);
        for (int j = 0; j < n; ++j) {
            if (xrow_[j] < 0 || xrow_[j] >= n || orow_[j] < 0 || orow_[j] >= n)
                throw std::invalid_argument("grid marker row out of range");
            if (xrow_[j] == orow_[j])
                throw std::invalid_argument("column has X and O in the same row");
            if (seenX[xrow_[j]]++ || seenO[orow_[j]]++)
                throw std::invalid_argument("grid markers do not form permutations");
        }
    }
};
