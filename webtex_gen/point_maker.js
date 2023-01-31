function point_maker() {
    this._point_helper = {
        _debug_tostr_from_ps: function (ps, fixed) {
            var outs = "";

            ps.forEach(function (v, idx) {
                if (idx != 0) outs += '\n';
                v.forEach(
                    function (p, idx) {
                        if (idx != 0) outs += ',';
                        outs += (
                            '(' + p.x.toFixed(fixed) + ',' +
                            p.y.toFixed(fixed) + ')'
                        );
                    }
                );
            });
            return outs;
        },

        /**
         * 
         * @param {Array} ps 需要疏散的点/点组 只能全为点或全为点组
         * @param {Number} mindx
         * @param {Number} mindy 
         * @returns 
         */
        _suppress_points: function (ps, mindx = - 1, mindy = -1) {
            var res = [];
            if (!Array.isArray(ps)) return [];
            if (ps.length < 3) return ps;
            if (Array.isArray(ps[0])) {
                ps.forEach(
                    function (v) {
                        res.push(_suppress_points(v, mindx, mindy));
                    }
                );
            }
            else ps.forEach(
                function (v, idx) {
                    if (idx == 0) {
                        res.push(v);
                        return;
                    }
                    var plast = res[res.length - 1];
                    if (Math.abs(v.x - plast.x) > mindx &&
                        Math.abs(v.y - plast.y) > mindy) res.push(v);
                }
            )
            return res;
        },
    };
    this._simplify_helper = {

        /**
         * Ramer-Douglas-Peucker 算法压缩/简化路径
         * @param {Array} points 
         * @param {Number} epsilon 
         * @returns 
         */
        _ramer_douglas_peucker: function (points, epsilon) {
            function _perpendicular_distance(point, start, end) {
                var area = Math.abs(0.5 * (start.x * end.y + end.x * point.y + point.x * start.y - end.x * start.y - point.x * end.y - start.x * point.y));
                var bottom = Math.sqrt(Math.pow(end.y - start.y, 2) + Math.pow(end.x - start.x, 2));
                var height = area / bottom * 2;
                return height;
            }
            function ramer_douglas_peucker(points, epsilon) {
                var first_point = points[0];
                var last_point = points[points.length - 1];
                if (points.length < 3) {
                    return points;
                }
                var index = -1;
                var dist = 0;
                for (var i = 1; i < points.length - 1; i++) {
                    var c_dist = _perpendicular_distance(points[i], first_point, last_point);
                    if (c_dist > dist) {
                        dist = c_dist;
                        index = i;
                    }
                }
                if (dist > epsilon) {
                    var l1 = points.slice(0, index + 1);
                    var l2 = points.slice(index);
                    var r1 = ramer_douglas_peucker(l1, epsilon);
                    var r2 = ramer_douglas_peucker(l2, epsilon);
                    var rs = r1.slice(0, r1.length - 1).concat(r2);
                    return rs;
                } else {
                    return [first_point, last_point];
                }
            };
            return ramer_douglas_peucker(points, epsilon)
        }

    }

    this.points = [];
    /**
    * 通过参数曲线取样
    * @param {object} func 参数函数 需包含 fx 与 fy 函数
    * @param {number[]} t_range
    * @param {number} t_count 取样点数(默认500)
    * @returns 由多个 连续点构成的数组 所构成的数组（内部通过识别 NaN 进行分割）
    */
    this.get_parametric_plot_points = function (
        func, t_range = [0, Math.PI * 2], t_count = 500
    ) {
        var res = [[]];
        var cur_tmp = res[res.length - 1];

        var L = Math.min(t_range[0], t_range[1]);
        var R = Math.max(t_range[0], t_range[1]);
        var p;

        var dt = Math.max((R - L) / t_count, 1e-9);

        var is_safe_num = ((x) => (
            Number.isFinite(x) &&
            Number.MIN_SAFE_INTEGER < x &&
            x < Number.MAX_SAFE_INTEGER
        ));

        var is_break_point = (x, f) => (
            Math.abs(f(x - Number.EPSILON) - f(x + Number.EPSILON)) >
            1e5 * Number.EPSILON
        );

        if (func.fx == undefined) {
            if (func.fy == undefined) return;
            func.fx = (t) => t;
        }
        if (func.fy == undefined) func.fy = (t) => t;

        while (L + dt < R) {
            if (is_break_point(L, func.fx) || is_break_point(L, func.fy)) {
                if (cur_tmp.length != 0) {
                    res.push([]);
                    cur_tmp = res[res.length - 1];
                }
                L += dt;
                continue;
            }
            p = { x: func.fx(L), y: func.fy(L) };
            L += dt;
            if (is_safe_num(p.x) && is_safe_num(p.y)) {
                cur_tmp.push(p);
            }
        }
        if (cur_tmp.length == 0) res.pop();
        return res;
    }

    // 通过极坐标下函数
    this.get_polar_plot_func = function (func, t_range, t_count) {
        if (func.theta == undefined) {
            if (func.rho == undefined) return [];
            func.theta = (t) => t;
        }
        if (func.rho == undefined) func.rho = (t) => t;
        return this.get_parametric_plot_points(
            {
                fx: (t) => func.rho(t) * Math.cos(func.theta(t)),
                fy: (t) => func.rho(t) * Math.sin(func.theta(t)),
            },
            t_range, t_count
        );
    }

    /**
     * 
     * @param {Object} func 包含 fx 或 fy 则当作参数方程，否则当作极坐标方程
     * @param {Array} t_range 范围 [Left,Right]
     * @param {Number} t_count 取采样点数
     */
    this.from = function (func, t_range, t_count) {
        if (typeof (func) == "function") {
            this.from({ fy: func }, t_range, t_count);
            return;
        }

        if ('fx' in func || 'fy' in func) {
            this.points = this.get_parametric_plot_points(func, t_range, t_count);
        }
        else
            this.points = this.get_polar_plot_func(func, t_range, t_count);
    }

    // simplify

    /**
     * 
     * @param {Number} epsilon 
     * @returns 
     */
    this.get_simplify_RDP = function (epsilon) {
        var tmp = [];
        var A = this._simplify_helper._ramer_douglas_peucker;
        this.points.forEach(function (v) {
            tmp.push(A(v, epsilon));
        })
        return tmp;
    }
    /**
     * 使用 Ramer-Douglas-Peucker 算法简化
     * @param {Number} epsilon 
     */
    this.simplify_RDP = function (epsilon) {
        this.points = this.get_simplify_RDP(epsilon);
    }

    // suppress

    this.get_suppress = function (mindx, mindy) {
        return this._point_helper._suppress_points(this.points, mindx, mindy);
    }
    /**
     * 疏散点
     * @param {number} mindx 最小x间隔
     * @param {number} mindy 最小y间隔
     */
    this.suppress = function (mindx, mindy) {
        this.points = this.get_suppress(mindx, mindy);
    }

    // toString

    this.toString = function (fixed = 3) {
        return this._point_helper._debug_tostr_from_ps(this.points, fixed);
    }
}

function catx(t) {
    function sin(x) { return Math.sin(x) };
    function cos(x) { return Math.cos(x) };
    return -((721 * sin(t)) / 4) + 196 / 3 * sin(2 * t) - 86 / 3 * sin(3 * t) - 131 / 2 * sin(4 * t) + 477 / 14 * sin(5 * t) + 27 * sin(6 * t) - 29 / 2 * sin(7 * t) + 68 / 5 * sin(8 * t) + 1 / 10 *
        sin(9 * t) + 23 / 4 * sin(10 * t) - 19 / 2 * sin(12 * t) - 85 / 21 * sin(13 * t) + 2 / 3 * sin(14 * t) + 27 / 5 * sin(15 * t) + 7 / 4 * sin(16 * t) + 17 / 9 * sin(17 * t) - 4 * sin(18 * t) -
        1 / 2 * sin(19 * t) + 1 / 6 * sin(20 * t) + 6 / 7 * sin(21 * t) - 1 / 8 * sin(22 * t) + 1 / 3 * sin(23 * t) + 3 / 2 * sin(24 * t) + 13 / 5 * sin(25 * t) + sin(26 * t) - 2 * sin(27 * t) +
        3 / 5 * sin(28 * t) - 1 / 5 * sin(29 * t) + 1 / 5 * sin(30 * t) + (2337 * cos(t)) / 8 - 43 / 5 * cos(2 * t) + 322 / 5 * cos(3 * t) - 117 / 5 * cos(4 * t) - 26 / 5 * cos(5 * t) - 23 / 3 *
        cos(6 * t) + 143 / 4 * cos(7 * t) - 11 / 4 * cos(8 * t) - 31 / 3 * cos(9 * t) - 13 / 4 * cos(10 * t) - 9 / 2 * cos(11 * t) + 41 / 20 * cos(12 * t) + 8 * cos(13 * t) + 2 / 3 * cos(14 * t) +
        6 * cos(15 * t) + 17 / 4 * cos(16 * t) - 3 / 2 * cos(17 * t) - 29 / 10 * cos(18 * t) + 11 / 6 * cos(19 * t) + 12 / 5 * cos(20 * t) + 3 / 2 * cos(21 * t) + 11 / 12 * cos(22 * t) - 4 / 5 *
        cos(23 * t) + cos(24 * t) + 17 / 8 * cos(25 * t) - 7 / 2 * cos(26 * t) - 5 / 6 * cos(27 * t) - 11 / 10 * cos(28 * t) + 1 / 2 * cos(29 * t) - 1 / 5 * cos(30 * t);
}
function caty(t) {
    function sin(x) { return Math.sin(x) };
    function cos(x) { return Math.cos(x) };
    return -((637 * sin(t)) / 2) - 188 / 5 * sin(2 * t) - 11 / 7 * sin(3 * t) - 12 / 5 * sin(4 * t) + 11 / 3 * sin(5 * t) - 37 / 4 * sin(6 * t) + 8 / 3 * sin(7 * t) + 65 / 6 * sin(8 * t) - 32 / 5 *
        sin(9 * t) - 41 / 4 * sin(10 * t) - 38 / 3 * sin(11 * t) - 47 / 8 * sin(12 * t) + 5 / 4 * sin(13 * t) - 41 / 7 * sin(14 * t) - 7 / 3 * sin(15 * t) - 13 / 7 * sin(16 * t) + 17 / 4 *
        sin(17 * t) - 9 / 4 * sin(18 * t) + 8 / 9 * sin(19 * t) + 3 / 5 * sin(20 * t) - 2 / 5 * sin(21 * t) + 4 / 3 * sin(22 * t) + 1 / 3 * sin(23 * t) + 3 / 5 * sin(24 * t) - 3 / 5 * sin(25 * t) +
        6 / 5 * sin(26 * t) - 1 / 5 * sin(27 * t) + 10 / 9 * sin(28 * t) + 1 / 3 * sin(29 * t) - 3 / 4 * sin(30 * t) - (125 * cos(t)) / 2 - 521 / 9 * cos(2 * t) - 359 / 3 * cos(3 * t) + 47 / 3 *
        cos(4 * t) - 33 / 2 * cos(5 * t) - 5 / 4 * cos(6 * t) + 31 / 8 * cos(7 * t) + 9 / 10 * cos(8 * t) - 119 / 4 * cos(9 * t) - 17 / 2 * cos(10 * t) + 22 / 3 * cos(11 * t) + 15 / 4 * cos(12 * t) -
        5 / 2 * cos(13 * t) + 19 / 6 * cos(14 * t) + 7 / 4 * cos(15 * t) + 31 / 4 * cos(16 * t) - cos(17 * t) + 11 / 10 * cos(18 * t) - 2 / 3 * cos(19 * t) + 13 / 3 * cos(20 * t) - 5 / 4 * cos(21 * t) +
        2 / 3 * cos(22 * t) + 1 / 4 * cos(23 * t) + 5 / 6 * cos(24 * t) + 3 / 4 * cos(26 * t) - 1 / 2 * cos(27 * t) - 1 / 10 * cos(28 * t) - 1 / 3 * cos(29 * t) - 1 / 19 * cos(30 * t);
}

function scan_line_fill(points, step, draw_line_func) {
    var nets = [], aets = [];
    var miny = Number.MAX_VALUE, maxy = -Number.MAX_VALUE;

    points.forEach(function (v) {
        if (v.y < miny) miny = v.y;
        if (v.y > maxy) maxy = v.y;
    });
    var line_id_of_y = (y) => Math.floor((y - miny) / step);
    var p_to_linep=(p)=>({y:line_id_of_y(p.y)*step+miny,x:p.x});
    var line_cnt = line_id_of_y(maxy) + 1;
    nets.length = line_cnt;

    //points.forEach((v)=>v.y=line_id_of_y(v.y)*step+miny);

    for (var i = 0; i < line_cnt; ++i) {
        nets[i] = [];
        points.forEach(function (cur_p, idx) {
            if(line_id_of_y(cur_p.y)!=i)return;
            var pre_p = p_to_linep(points[(idx - 1 + points.length) % points.length]);
            var nxt_p = p_to_linep(points[(idx + 1 + points.length) % points.length]);

            if (pre_p.y > cur_p.y) {
                nets[i].push({
                    x: cur_p.x,
                    dx: (pre_p.x - cur_p.x) / (pre_p.y - cur_p.y),
                    ymax: pre_p.y
                });
            }
            if (nxt_p.y > cur_p.y) {
                nets[i].push({
                    x: cur_p.x,
                    dx: (nxt_p.x - cur_p.x) / (nxt_p.y - cur_p.y),
                    ymax: nxt_p.y
                });
            }
        });
    }

    for (var i = 0; i < line_cnt; ++i) {
        aets.forEach((v) => (v.x += v.dx * step));
        aets = aets.filter((v) => (line_id_of_y(v.ymax)>i));

        nets[i].forEach((v) => (aets.push(v)));
        aets.sort((a, b) => (a.x - b.x));

        for (var j = 0; j + 1 < aets.length; j += 2) {
            draw_line_func(aets[j].x, i * step + miny, aets[j + 1].x, i * step + miny + step);
        }
    }
}

var ps_mker = new point_maker();
//console.log(ps_mker.toString());
function F(fixed, step,zoom) {
    var outs = "";
    var cnt = 0;
    ps_mker.points.forEach((v) => {
        v.forEach((p) => { p.x *= zoom, p.y *= zoom; });
        scan_line_fill(v, step, (sx, sy, ex, ey) => {
            var _my_fixed=(n)=>parseFloat(n.toFixed(fixed));
            cnt++;
            outs += "\\kern" + _my_fixed(sx) + "pt\\rule[" + _my_fixed(sy) + "pt]" 
                + _my_fixed(ex - sx) + "pt"
                + _my_fixed(ey - sy) + "pt\\kern"
                + _my_fixed(-ex) + "pt";
        });
        v.forEach((p) => { p.x /= zoom, p.y /= zoom; });
    })
    console.log("fill line count:",cnt);
    return outs;
}
