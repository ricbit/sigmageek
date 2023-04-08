#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <queue>
#include <array>
#include <functional>
#include <map>
#include <numeric>
#include <algorithm>
#include <execution>
#include <fstream>

using namespace std;
using namespace std::execution;
using uint8 = unsigned char;

struct Position {
    int y : 12;
    int x : 13;
    unsigned int direction : 2;
    unsigned int level: 16;
    unsigned int prev;
    unsigned int children: 3;

    bool same(const Position& other) const {
        return this->x == other.x and this->y == other.y;
    }
};

constexpr int msize = 2048;
constexpr int score_num = 8, score_den = 5;

struct Grid {
    Grid(const vector<bool>& grid, int sy, int sx, Position start, Position end)
        : grid(grid), sy(sy), sx(sx), start(start), end(end) {
    }

    void print() {
        for (int j = 0; j < sy; j++) {
            for (int i = 0; i < sx; i++) {
                cout << get(j, i);
                if (i < sx - 1) {
                    cout << " ";
                }
            }
            cout << "\n";
        }
        cout << "\n";
    }

    const int get(int y, int x) const {
        return grid[y * msize + x] ? 1 : 0;
    }

    const vector<bool> grid;
    const int sy, sx;
    const Position start, end;
};

void trim_back(string& line) {
    while (line.back() <= 32) {
        line.pop_back();
    }
}

Grid read_input() {
    vector<bool> temp_grid(msize * msize, false);
    string line;
    int rows = 0, cols = 0;
    Position start{0, 0, 0, 0, 0, 0}, end{0, 0, 0, 0, 0, 0};
    for (;getline(cin, line); rows++) {
        trim_back(line);
        vector<uint8> row;
        row.reserve(msize);
        stringstream iss(line);
        int value;
        while (!iss.eof()) {
            iss >> value;
            row.push_back(value);
        }
        cols = row.size();
        for (int i = 0; i < cols; i++) {
            temp_grid[rows * msize + i] = row[i] == 1;
            if (row[i] == 3) {
                start.y = rows;
                start.x = i;
            }
            if (row[i] == 4) {
                end.y = rows;
                end.x = i;
            }
        }
    }
    Grid grid(temp_grid, rows, cols, start, end);
    return grid;
}

vector<char> read_movements(string filename) {
    ifstream f(filename, ios_base::in);
    vector<char> movs;
    char c;
    f >> c;
    while (!f.eof()) {
        char d;
        f >> d;
        movs.push_back(c);
        c = d;
    }
    cerr << movs.size() << "\n";
    return movs;
}

bool valid_coord(const Grid& grid, int y, int x) {
    return x >= 0 and x < grid.sx and y >= 0 and y < grid.sy; 
}

bool valid_slow(const Grid& grid, int jj, int ii) {
    return valid_coord(grid, jj, ii) and grid.get(jj,ii) == 1;
}

bool valid_fast(const Grid& grid, int jj, int ii) {
    return grid.get(jj,ii) == 1;
}

template<typename T, typename V>
void sum_neigh(const Grid& grid, int y, int x, T action, V validation) {
    for(int j = -1; j <= 1; j++) {
        for (int i = -1; i <= 1; i++) {
            if (i == 0 and j == 0) {
                continue;
            }
            int jj = y + j;
            int ii = x + i;
            if (validation(grid, jj, ii)) {
                action(grid.get(jj,ii));
            }
        }
    }
}

uint8 transition_rule(const Grid& grid, int y, int x, int block) {
    if (grid.get(y,x) == 0 and block > 1 and block < 5) {
        return 1;
    } else if (grid.get(y,x) == 1 and (block <= 3 or block >= 6)) {
        return 0;
    }
    return grid.get(y,x);
}

template<typename V>
uint8 transition(const Grid &grid, int j, int i, V validation) {
    int block = 0;
    sum_neigh(grid, j, i, [&block](int x){ block += x; }, validation);
    return transition_rule(grid, j, i, block);
}

template<typename T>
void iter_line(vector<bool>& newgrid, int j, int i1, int i2, T callback) {
    for (int i = i1; i < i2; i++) {
        newgrid[j * msize + i] = callback(j, i);
    }
}

Grid iter_grid(const Grid& grid, const vector<int>& yrange) {
    vector<bool> newgrid(grid.grid);
    auto slow = [&grid](int j, int i) { return transition(grid, j, i, valid_slow); };
    auto fast = [&grid](int j, int i) { return transition(grid, j, i, valid_fast); };
    iter_line(newgrid, 0, 0, grid.sx, slow);
    for_each(std::execution::par, yrange.begin(), yrange.end(), [&grid, &newgrid, &slow, &fast](int j){
        iter_line(newgrid, j, 0, 1, slow);
        iter_line(newgrid, j, 1, grid.sx - 1, fast);
        iter_line(newgrid, j, grid.sx - 1, grid.sx, slow);
    });
    iter_line(newgrid, grid.sy - 1, 0, grid.sx, slow);
    newgrid[grid.start.y * msize + grid.start.x] = false;
    newgrid[grid.end.y * msize + grid.end.x] = false;
    return Grid(newgrid, grid.sy, grid.sx, grid.start, grid.end);
}

struct Direction {
    int y, x;
    char name;
};

constexpr array<Direction, 4> directions = {{
    {-1,  0, 'U'},
    { 1,  0, 'D'}, 
    { 0, -1, 'L'}, 
    { 0,  1, 'R'}
}};

struct Head {
    int score;
    unsigned int position;
    auto operator<(const Head& other) const {
        return tie(score, position) < tie(other.score, other.position);
    }
};

int score(const Grid& grid, const Position& pos, const Position &goal) {
    return -(pos.level + score_num * (abs(goal.x - pos.x) + abs(goal.y - pos.y)) / score_den);
}

string path(const vector<Position>& positions, unsigned int start) {
    string result;
    for (auto cur = start; !(cur == 0 && positions[cur].prev == 0); cur = positions[cur].prev) {
        result = directions[positions[cur].direction].name + " "s + result;
    }
    trim_back(result);
    return result;
}

const Grid& get_grid(vector<Grid>& grid_cache, unsigned int level, const vector<int>& yrange) {
    while (grid_cache.size() < level + 1) {
        grid_cache.push_back(iter_grid(*grid_cache.rbegin(), yrange));
    }
    return grid_cache[level];
}

struct BackTrie {
    BackTrie(Position start) : positions{start}, left(0), next(0) {}

    void pop(int pos) {
        left++;
        int parent = positions[pos].prev;
        positions[parent].children--;
        positions[pos].prev = next;
        next = pos;
        if (positions[parent].children == 0) {
            pop(parent);
        }
    }

    unsigned int push(Position value) {
        if (left) {
            left--;
            int cur = next;
            next = positions[cur].prev;
            positions[cur] = value;
            return cur;
        } else {
            positions.push_back(value);
            return positions.size() - 1;
        }
    }

    vector<Position> positions;
    int left, next;
};

struct VisitedPosition {
    int y : 12;
    int x : 12;
    int level : 16;

    auto operator<(const VisitedPosition& other) const {
        return tie(y, x, level) < tie(other.y, other.x, other.level);
    }
};

struct VisitedQuads {
    VisitedQuads() {}
    
    void set(int j, int i, int layer) {
        auto key = VisitedPosition{j >> 6, i >> 6, layer};
        auto it = visited.find(key);
        if (it == visited.end()) {
            auto it_pair = visited.emplace(key, vector<bool>(64 * 64, false));
            it = it_pair.first;
        }
        it->second[(j & 0x3F) * 64 + (i & 0x3F)] = true;
    }

    bool check(int j, int i, int layer) {
        auto key = VisitedPosition{j >> 6, i >> 6, layer};
        auto it = visited.find(key);
        if (it == visited.end()) {
            return false;
        }
        return it->second[(j & 0x3F) * 64 + (i & 0x3F)];
    }

    map<VisitedPosition, vector<bool>> visited;
};
string search(const Grid& grid) {
    BackTrie trie{grid.start};
    priority_queue<Head> next;

    vector<Grid> grid_cache{grid};
    VisitedQuads visited;
    visited.set(grid.start.y, grid.start.x, 0);
    next.push({score(grid, grid.start, grid.end), 0});
    vector<int> yrange(grid.sy - 2);
    iota(yrange.begin(), yrange.end(), 1);
    int steps = 0;
    while (!next.empty()) {
        Head current = next.top();
        next.pop();
        auto pos = trie.positions[current.position];
        if (pos.same(grid.end)) {
            cerr << pos.level << endl;
            return path(trie.positions, current.position);
        }
        if (++steps % 100 == 0) {
            cerr << steps << " " << trie.positions.size() << " " 
                 << -current.score << " " << pos.level << "\n";
        }
        unsigned int next_level = pos.level + 1;
        const Grid& next_grid = get_grid(grid_cache, next_level, yrange);
        trie.positions[current.position].children = 0;
        for (unsigned int d = 0; d < 4; d++) {
            int jj = pos.y + directions[d].y;
            int ii = pos.x + directions[d].x;
            if (valid_coord(grid, jj, ii)) {
                if (visited.check(jj, ii, next_level) or next_grid.get(jj,ii) == 1) {
                    continue;
                }
                Position next_position{jj, ii, d, next_level, current.position, 0};
                auto next_index = trie.push(next_position);
                visited.set(jj, ii, next_level);
                auto next_score = score(next_grid, next_position, grid.end);
                next.push({next_score, next_index});
                trie.positions[current.position].children++;
            }
        }
        if (trie.positions[current.position].children == 0) {
            trie.pop(current.position);
        }
    }
    return "nope\n"s;
}

const Direction& get_dir(char mov) {
    for (auto& dir : directions) {
        if (mov == dir.name) {
            return dir;
        }
    }
    cerr << "Invalid direction" << endl;
    exit(1);
}

void check(const Grid& grid, const vector<char>& movs) {
    vector<Grid> grid_cache{grid};
    Position pos = grid.start;
    vector<int> yrange(grid.sy - 2);
    iota(yrange.begin(), yrange.end(), 1);
    int lost = 0;
    for (int i = 0; i < static_cast<int>(movs.size()); i++) {
        if (i % 100 == 0) {
            cerr << i << endl;
        }
        auto current = get_grid(grid_cache, i + 1, yrange);
        auto dir = get_dir(movs[i]);
        pos = Position{pos.y + dir.y, pos.x + dir.x};
        if (!valid_coord(current, pos.y, pos.x)) {
            cerr << "Out of grid" << endl;
            exit(1);
        }
        if (current.get(pos.y, pos.x) == 1) {
            lost++;
            cerr << i << "\n";
            for (int jj = -2; jj <= 2; jj++) {
                for (int ii = -2; ii <= 2; ii++) {
                    cerr << get_grid(grid_cache, i, yrange).get(pos.y + jj, pos.x + ii) << " ";
                }
                cerr << "\n";
            }
            cerr << " to \n";
            for (int jj = -2; jj <= 2; jj++) {
                for (int ii = -2; ii <= 2; ii++) {
                    cerr << get_grid(grid_cache, i + 1, yrange).get(pos.y + jj, pos.x + ii) << " ";
                }
                cerr << "\n";
            }
            vector<bool> scratch = get_grid(grid_cache, i, yrange).grid;
            scratch[pos.y * msize + pos.x] = false;
            Grid newgrid = iter_grid(Grid(scratch, grid.sy, grid.sx, grid.start, grid.end), yrange);
            cerr << " fix \n";
            for (int jj = -2; jj <= 2; jj++) {
                for (int ii = -2; ii <= 2; ii++) {
                    cerr << newgrid.get(pos.y + jj, pos.x + ii) << " ";
                }
                cerr << "\n";
            }
        }
    }
    cerr << "\n" << movs.size() << " movements\n";
    cerr << lost << " lives lost\n";
    cerr << (grid.end.same(pos) ? "Succeeded\n" : "Failed\n");
}

int main(int argc, char **argv) {
    Grid grid = read_input();
    auto movs = read_movements(argv[1]);
    check(grid, movs);
    /*for (char c : movs) {
        cerr << "-" << c << "- ";
    }*/
    cerr << endl;
    return 0;
}