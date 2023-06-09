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
#include <set>

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

    bool operator<(const Position& other) const {
        return tie(y, x, level) < tie(other.y, other.x, other.level);
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

pair<vector<vector<char>>, vector<int>> read_movements(string filename) {
    ifstream f(filename, ios_base::in);
    vector<vector<char>> movs;
    vector<int> start;
    string line;
    for (;getline(f, line);) {
        trim_back(line);
        stringstream iss(line);
        int lag;
        iss >> lag;
        cerr << lag << endl;
        start.push_back(lag);
        vector<char> row;
        char c;
        iss >> c;
        while (!iss.eof()) {
            char d;
            iss >> d;
            row.push_back(c);
            c = d;
        }
        movs.push_back(row);
    }
    cerr << movs.size() << "\n";
    return make_pair(movs, start);
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

const Grid& get_grid(vector<Grid>& grid_cache, unsigned int level, const vector<int>& yrange) {
    while (grid_cache.size() < level + 1) {
        grid_cache.push_back(iter_grid(*grid_cache.rbegin(), yrange));
    }
    return grid_cache[level];
}

struct VisitedPosition {
    int y : 12;
    int x : 12;
    int level : 16;

    auto operator<(const VisitedPosition& other) const {
        return tie(y, x, level) < tie(other.y, other.x, other.level);
    }
};

const Direction& get_dir(char mov) {
    for (auto& dir : directions) {
        if (mov == dir.name) {
            return dir;
        }
    }
    cerr << "Invalid direction" << endl;
    exit(1);
}

void check(const Grid& grid, const vector<vector<Position>>& particles) {
    vector<Grid> grid_cache{grid};
    Position pos = grid.start;
    vector<int> yrange(grid.sy - 2);
    iota(yrange.begin(), yrange.end(), 1);
    for (int p = 0; p < static_cast<int>(particles.size()); p++) {
        int lost = 0;
        for (int i = 0; i < static_cast<int>(particles[p].size()); i++) {
            pos = particles[p][i];
            auto current = get_grid(grid_cache, pos.level + 1, yrange);
            if (!valid_coord(current, pos.y, pos.x)) {
                cerr << "Out of grid" << endl;
                exit(1);
            }
            if (current.get(pos.y, pos.x) == 1) {
                lost++;
            }
        }
        cerr << "\n" << particles[p][0].level << " " << particles[p].size() << " movements\n";
        cerr << lost << " lives lost\n";
        cerr << (grid.end.same(pos) ? "Succeeded\n" : "Failed\n");
    }
}

void collision(const vector<vector<Position>>& particles) {
    set<Position> unique;
    int total = 0;
    for (auto& particle : particles) {
        for (auto& pos : particle) {
            Position scrub = pos;
            scrub.children = 0;
            scrub.direction = 0;
            scrub.prev = 0;
            unique.insert(scrub);
            total++;
        }
    }
    cerr << "Collisions: " << (total - static_cast<int>(unique.size())) << "\n";
}

vector<vector<Position>> gen_positions(
        const vector<vector<char>>& movs, const vector<int>& lags, Position start) {
    vector<vector<Position>> positions;
    for (int lag = 0; lag < static_cast<int>(movs.size()); lag++) {
        Position pos = start;
        vector<Position> path;
        for (int p = 0; p < static_cast<int>(movs[lag].size()); p++) {
            pos.y += get_dir(movs[lag][p]).y;
            pos.x += get_dir(movs[lag][p]).x;
            pos.level = lags[lag] + p;
            path.push_back(pos);
        }
        positions.push_back(path);
    }
    return positions;
}

int main(int argc, char **argv) {
    Grid grid = read_input();
    auto pmovs = read_movements(argv[1]);
    auto& movs = pmovs.first;
    auto& lags = pmovs.second;
    for (int i = 0; i < static_cast<int>(movs.size()); i++) {
        cerr << lags[i] << " " << movs[i].size() << "\n";
    }
    auto particles = gen_positions(movs, lags, grid.start);
    check(grid, particles);
    collision(particles);
    cerr << "Total particles:  " << movs.size() << endl;
    return 0;
}