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
#include <unordered_map>

using namespace std;
using namespace std::execution;
using uint8 = unsigned char;


constexpr int msize = 512;
constexpr int score_num = 1, score_den = 1;
constexpr int initial_lives = 0;
constexpr int size_bits = 10;

struct Position {
    int y : size_bits;
    int x : size_bits;
    unsigned int direction : 2;
    unsigned int level : 16;
    unsigned int prev;
    unsigned int children : 3;

    bool same(const Position& other) const {
        return this->x == other.x and this->y == other.y;
    }
};

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
    for_each(std::execution::par_unseq, yrange.begin(), yrange.end(), [&grid, &newgrid, &slow, &fast](int j){
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

int score(const Grid& grid, const Position& pos, Position goal) {
    return -(pos.level + score_num * (abs(goal.x - pos.x) + abs(goal.y - pos.y)) / score_den);
}

vector<uint8> path(const vector<Position>& positions, unsigned int start) {
    vector<uint8> result;
    for (auto cur = start; !(cur == 0 && positions[cur].prev == 0); cur = positions[cur].prev) {
        result.push_back(positions[cur].direction);
    }
    reverse(result.begin(), result.end());
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
    int y : size_bits;
    int x : size_bits;
    unsigned int level : 16;

    auto operator<(const VisitedPosition& other) const {
        return tie(y, x, level) < tie(other.y, other.x, other.level);
    }

    auto operator==(const VisitedPosition& other) const {
        return tie(y, x, level) == tie(other.y, other.x, other.level);
    }
};

struct VisitedPositionHash {
    auto operator()(const VisitedPosition& pos) const {
        return std::hash<int>()(pos.y) ^ std::hash<int>()(pos.x) ^ std::hash<unsigned int>()(pos.level);
    }
};

struct VisitedQuads {
    VisitedQuads() {}
    
    void set(int j, int i, unsigned int layer) {
        auto key = VisitedPosition{j >> 6, i >> 6, layer};
        auto it = visited.find(key);
        if (it == visited.end()) {
            auto it_pair = visited.emplace(key, vector<bool>(64 * 64, false));
            it = it_pair.first;
        }
        it->second[(j & 0x3F) * 64 + (i & 0x3F)] = true;
    }

    bool check(int j, int i, unsigned int layer) {
        auto key = VisitedPosition{j >> 6, i >> 6, layer};
        auto it = visited.find(key);
        if (it == visited.end()) {
            return false;
        }
        return it->second[(j & 0x3F) * 64 + (i & 0x3F)];
    }

    unordered_map<VisitedPosition, vector<bool>, VisitedPositionHash> visited;
};

int distance(const Position& pos, const Position& goal) {
    return abs(goal.x - pos.x) + abs(goal.y - pos.y);
}

bool block(const Position& pos, const Position& goal, int lag) {
    return pos.level < 49900 and pos.same(goal);
}

void trim_heap(priority_queue<Head>& next) {
    unsigned int prune = 100;
    unsigned int keep = 40;
    if (next.size() > prune) {
        vector<Head> save;            
        while (next.size() > prune - keep) {
            save.push_back(next.top());
            next.pop();
        }
        while (!next.empty()) {
            next.pop();
        }
        for (auto& head : save) {
            next.push(head);
        }
    }
}

vector<uint8> search(const Grid& grid, vector<Grid>& grid_cache, VisitedQuads& global, int lag) {
    VisitedQuads visited;
    Position lag_start = grid.start;
    lag_start.level = lag;
    BackTrie trie{lag_start};
    priority_queue<Head> next;
    visited.set(lag_start.y, lag_start.x, 0);
    next.push({score(grid, lag_start, grid.end), 0});
    vector<int> yrange(grid.sy - 2);
    iota(yrange.begin(), yrange.end(), 1);
    int steps = 0;
    while (!next.empty()) {
        Head current = next.top();
        next.pop();
        trim_heap(next);
        auto pos = trie.positions[current.position];
        if (pos.same(grid.end)) {
            cerr << pos.level << endl;
            return path(trie.positions, current.position);
        }
        if (++steps % 100000 == 0) {
            cerr << steps << " " << trie.positions.size() << " " 
                 << -current.score << " " << pos.level << endl;
        }
        unsigned int next_level = pos.level + 1;
        const Grid& next_grid = get_grid(grid_cache, next_level, yrange);
        trie.positions[current.position].children = 0;
        for (unsigned int d = 0; d < 4; d++) {
            int jj = pos.y + directions[d].y;
            int ii = pos.x + directions[d].x;
            if (valid_coord(grid, jj, ii) and !block(Position{jj, ii, d, next_level}, grid.end, lag)) {
                if (visited.check(jj, ii, next_level) 
                    or global.check(jj, ii, next_level)
                    or next_grid.get(jj, ii) == 1) {
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
    return vector<uint8>();
}

int main() {
    Grid grid = read_input();
    vector<Grid> grid_cache{grid};
    vector<vector<uint8>> particles;
    int total = 0;
    VisitedQuads global;
    for (int i = 0; i < 40000; i++) {
        particles.emplace_back(search(grid, grid_cache, global, i));
        auto& current = *particles.rbegin();
        if (!current.empty() && current.size() < 49999) {
            cerr << "total " << ++total << "\n";
            cout << i << " ";
            Position pos = grid.start;
            for (int p = 0; p < static_cast<int>(current.size()); p++) {
                global.set(pos.y, pos.x, p + i);
                pos.y += directions[current[p]].y;
                pos.x += directions[current[p]].x;
                cout << directions[current[p]].name;
                if (p < static_cast<int>(current.size()) - 1) {
                    cout << " ";
                }
            }
            cout << endl;
        }
    }
    return 0;
}