// you can test individual functions in the Console/Web developer tools in either Firefox or Chrome
// the "number" type is a double


// Ensure that KaTeX is loaded before running this.
declare var katex: any;

// Assert that a condition is true, otherwise crash the program with a message.
function assert(cond: boolean, msg?: string) {
    if (cond)
        return;

    const message = "Assertion failed" + ((msg == undefined) ? "." : ": " + msg);
    throw new Error(message);
}

// Memoise a function
//
// Assume that the given function takes values in "K" and returns a value of type "V"
// this caches the result of calls so that they do not need to be computed twice
//
// You do not really need it; create a class SymmetricIrrepInfo that takes a single parameter "n"
// that will store the information for Sn. Then, at construction time, compute the needed factorials,
// conjugacy classes, etc...
function memoised<K, V>(f: (k: K) => V): (k: K) => V {
    let memoTable: {[s: string]: V} = {};
    return function(k: K): V {
        let kString = JSON.stringify(k);
        if (!memoTable.hasOwnProperty(kString))
            memoTable[kString] = f(k);
        return memoTable[kString];
    }
}

// Generate the list of numbers [1, 2, ..., n].
function range(n: number): number[] {
    let numbers: number[] = [];
    for (let i = 1; i <= n; i++)
        numbers.push(i);

    return numbers;
}

// The sum of a list.
function sum(ns: number[]): number {
    return ns.reduce((a, b) => a + b, 0);
}

// The product of a list.
function product(ns: number[]): number {
    return ns.reduce((a, b) => a * b, 1);
}

// The factorial of n.
let factorial = memoised(factorialHelp);
function factorialHelp(n: number): number {
    assert(n >= 0);
    return product(range(n));
}

// The greatest common divisor of two integers.
function gcd(a: number, b: number): number {
    if (b == 0)
        return a;
    return gcd(b, a % b);
}


// ------------------
// --- Partitions ---
// ------------------

// A partition is a weakly decreasing list of positive integers.
type Partition = number[]; // so... a double row vector

// Check that a list is weakly decreasing, i.e. [a >= b >= ... >= d].
function isWeaklyDecreasing(ns: number[]): boolean {
    for (let i = 0; i < ns.length - 1; i++)
        if (ns[i] < ns[i + 1])
            return false;

    return true;
}

// Check a list is a partition.
function isPartition(ns: number[]): boolean {
    return isWeaklyDecreasing(ns) && ns.every(x => x > 0);
}

// Turn any weak composition (list of nonnegative integers) into a partition,
// by removing any zero parts, and sorting.
function partitionFromComposition(comp: number[]): Partition {
    return comp.filter(x => x > 0).sort((a, b) => b - a);
}

// Generate all partitions of n, in lex order.
let partitionsOf = memoised(partitionsOfHelper);
function partitionsOfHelper(n: number): Partition[] {
    // Helper function: returns all the partitions of n which fit inside a width w strip.
    function helper(n: number, w: number): Partition[] {
        assert(n >= 0);
        assert(w >= 0);

        if (n == 0)
            return [[]];
        if (w <= 0)
            return [];

        let results: Partition[] = [];
        for (let row = Math.min(n, w); row >= 1; row--) {
            let smaller = helper(n - row, row);
            results.push(...smaller.map(part => [row].concat(part)));
        }

        return results;
    }
    return helper(n, n);
}

// Turn a partition into "group notation", for example the partition
// [4, 4, 3, 2, 2, 2, 1] => [[4, 2], [3, 1], [2, 3], [1, 1]].
function groupPartition(partition: Partition): [number, number][] {
    if (partition.length == 0)
        return [];

    let i = 0, j = 0;
    let groups: [number, number][] = [];
    while (i < partition.length) {
        for (j = i; j < partition.length; j++)
            if (partition[i] != partition[j])
                break;

        groups.push([partition[i], j - i]);
        i = j;
    }

    return groups;
}

// The symmetric group S_l acts on a partition of length l by permuting the parts.
// When the partition is written in "group notation" [[a1, b1], ..., [ar, br]], the
// order of the stabiliser of this action is b1! * ... * br!.
function stabiliserOrder(partition: Partition): number {
    // MATLAB: .map is the equivalent of arrayfun or cellfun
    let factorials = groupPartition(partition).map(([_, b]) => factorial(b));
    return product(factorials);
}

// For a partition of n, how many elements of the symmetric group S_n are conjugate
// to an element with the given cycle type? When the cycle type is written in
// "group notation" [[a1, b1], ..., [ar, br]], the answer is
// n! / ((a1^b1 * ... * ar^br) * (b1! * ... * br!)
let conjugacySize = memoised(conjugacySizeHelp);
function conjugacySizeHelp(cycleType: Partition): number {
    // MATLAB: just use the standard formula, this is what they do!
    let n = sum(cycleType);
    let groups = groupPartition(cycleType); // MATLAB: this is a list of tuples
    let factorials = groups.map(([_, b]) => factorial(b)); // the map function maps every tuple to something, here only the exponents are kept
    let powers = groups.map(([a, b]) => Math.pow(a, b)); // same here

    return factorial(n) / (product(powers) * product(factorials));
}

// For an element of the symmetric group with a given cycle type, what cycle type
// does the pth power of that element have? The answer can be computed independently
// for each cycle; for a cycle of size n, the power n^p breaks into gcd(n, p) cycles,
// each of length n / gcd(n, p).
function cyclePow(cycleType: Partition, p: number): Partition {
    let newCycleType: number[] = [];
    for (let cycle of cycleType) {
        let div = gcd(cycle, p);
        for (let i = 0; i < div; i++)
            newCycleType.push(cycle / div); // MATLAB: newCycleType(end+1) = cycle/div
    }
    return partitionFromComposition(newCycleType);
}

// Turn a partition into a string
function showPartition(partition: Partition): string {
    return JSON.stringify(partition); // MATLAB: here you need to be more clever:  ['p_' strjoin(arrayfun(@num2str, partition, 'uniform', 0), '_')]
}

// Turn a string into a partition
function readPartition(partition: string): Partition {
    return JSON.parse(partition); // res = cellfun(@str2num, strsplit(partition(3:end), '_'))
}


// MATLAB: we don't care of this one

// Return a LaTeX expression of the partition in "group notation", so
// [4, 4, 3, 2, 2, 2, 1] => [4^2, 3^1, 2^3, 1^1].
function showPartitionCompact(partition: Partition): string {
    return "[" + groupPartition(partition).map(([a, b]) => `${a}^{${b}}`).join(",") + "]";
}


// MATLAB: just write a class using a "store" property which is a struct
// MATLAB: don't care about element types, it's never a thing, MATLAB is dynamically typed

// A map of type (Partition -> T). Since Javascript only understands strings as
// may keys, we convert to and from strings a lot.
class PartMap<T> {
    private store: {[s: string]: T} = {} // MATLAB: this is going to be a struct; index with "store.(fieldName)"; check presernve "isfield(store, fieldName)"

    has(k: Partition): boolean {
        return this.store.hasOwnProperty(showPartition(k));
    }

    get(k: Partition): (null | T) {
        let kString = showPartition(k);
        if (!this.store.hasOwnProperty(kString))
            return null;

        return this.store[kString];
    }

    set(k: Partition, v: T): void {
        this.store[showPartition(k)] = v;
    }

    keys(): Partition[] {
        return Object.keys(this.store).map(readPartition);
    }
}

// -----------------------------------------
// --- Linear combinations of partitions ---
// -----------------------------------------

// MATLAB: same here, just inherit the PartMap class

// An integer-linear combination of partitions.
class Lin extends PartMap<number> {
    get(k: Partition): number {
        let coeff = super.get(k);
        if (coeff == null)
            return 0;
        return coeff;
    }

    support(): Partition[] {
        return this.keys().filter(k => this.get(k) != 0);
    }

    toList(): [Partition, number][] {
        return this.support().map(k => [k, this.get(k)]);
    }

    show(letter?: string): string {
        const output: string[] = [];
        for (const part of this.support()) {
            let coeffStr = "" + this.get(part);
            let partStr = showPartition(part);

            if (partStr == "[]") {
                output.push(coeffStr);
                continue;
            }
            if (coeffStr == "1")
                coeffStr = "";

            output.push(coeffStr + ((letter === undefined) ? '' : letter) + partStr)
        }
        if (output.length == 0)
            return "0";
        return output.join(" + ");
    }

    static fromList(pairs: [Partition, number][]): Lin {
        let lin = new Lin();
        for (let [partition, coeff] of pairs)
            lin.set(partition, lin.get(partition) + coeff);

        return lin;
    }

    static zero(): Lin {
        return Lin.fromList([]);
    }

    static scale(scalar: number, lin: Lin): Lin {
        return Lin.fromList(lin.toList().map(([k, v]) => [k, scalar * v]));
    }

    static add(left: Lin, right: Lin): Lin {
        return Lin.fromList(left.toList().concat(right.toList()));
    }

    static sub(left: Lin, right: Lin): Lin {
        return Lin.add(left, Lin.scale(-1, right));
    }

    static applyDiagonal(f: (a: Partition) => number, lin: Lin): Lin {
        return Lin.fromList(lin.toList().map(([k, v]) => [k, f(k) * v]));
    }

    static applyBilinear(f: (a: Partition, b: Partition) => Lin, left: Lin, right: Lin): Lin {
        let pairs: [Partition, number][] = [];
        for (let [k1, v1] of left.toList()) // MATLAB: same thing here as before, "left.toList" returns a list of pairs, and [k1, v1] destructures the pairs into two elements
            for (let [k2, v2] of right.toList())
                pairs.push(...Lin.scale(v1 * v2, f(k1, k2)).toList()); // this is just a concatenation of lists I guess

        return Lin.fromList(pairs);
    }
}

// ------------------------------------
// --- Power sum to monomial basis  ---
// ------------------------------------

// Multiply the monomial symmetric function m_(r) with an augmented monomial symmetric function.
// The product M_r * M(l1, ..., lk) is equal to the sum
//    M(l1 + r, ..., lk) + ... + M(l1, ..., lk + r) + M(r, l1, ..., lk),
// where we need to remember to sort the compositions on the right side back into being partitions.
function multRowWithAugmented(rowPartition: Partition, partition: Partition): Lin {
    if (partition.length == 1)
        [rowPartition, partition] = [partition, rowPartition];
    assert(rowPartition.length == 1);
    let row = rowPartition[0];
    let compositions = partition.map(
        (part, i) => partition.slice(0, i).concat([part + row]).concat(partition.slice(i + 1)));
    compositions.push(partition.concat([row]));
    return Lin.fromList(compositions.map(comp => [partitionFromComposition(comp), 1]));
}

// Expand a power sum symmetric function into an augmented monomial symmetric function.
function powerToAugmentedMonomial(partition: Partition): Lin {
    return partition
        .map(part => Lin.fromList([[[part], 1]]))
        .reduce((a, b) => Lin.applyBilinear(multRowWithAugmented, a, b), Lin.fromList([[[], 1]]));
}

// Expand a power sum symmetric function into monomial symmetric function.
function powerToMonomial(partition: Partition): Lin {
    return Lin.applyDiagonal(stabiliserOrder, powerToAugmentedMonomial(partition));
}

// ----------------------------------------
// --- Characters and character tables  ---
// ----------------------------------------

// A character of the symmetric group is a linear combination of partitions. A character table
// is a map from partitions to such linear combinations.

// MATLAB: This thing is fun: it is a map from partitions to (map from partitions to double)
class CharacterTable extends PartMap<Lin> {
    readonly partitions: Partition[]

    constructor(readonly n: number) {
        super();
        this.partitions = partitionsOf(n);
        for (let partition of this.partitions)
            this.set(partition, Lin.zero());
    }

    get(partition: Partition): Lin {
        let result = super.get(partition);
        assert(result != null);
        return <Lin>result;
    }
}

// To generate the character table for the permutation modules, we just need to transpose the
// turn-power-sums-into-monomial data.
function characterTablePermutation(n: number): CharacterTable {
    let table = new CharacterTable(n);

    for (let col of table.partitions)
        for (let [row, coeff] of powerToMonomial(col).toList())
            table.set(row, Lin.add(table.get(row), Lin.fromList([[col, coeff]])))

    return table;
}

// Inner product on class functions in the symmetric group on n letters.
function innerProduct(n: number, left: Lin, right: Lin): number {
    return sum(partitionsOf(n).map(part => left.get(part) * right.get(part) * conjugacySize(part))) / factorial(n);
}

// Character table for the Specht modules (we will memoise this).
function characterTableSpechtHelper(n: number): CharacterTable {
    let permTable = characterTablePermutation(n);
    let partitions = partitionsOf(n);
    for (let i = 0; i < partitions.length; i++) {
        // This row of permTable is currently a Specht module. We look at everything below it, compute
        // the inner product, and subtract off where needed.
        let spechtMod = permTable.get(partitions[i]);
        for (let j = i + 1; j < partitions.length; j++) {
            let mod = permTable.get(partitions[j]);
            let multiplicity = innerProduct(n, spechtMod, mod);
            permTable.set(partitions[j], Lin.sub(mod, Lin.scale(multiplicity, spechtMod)));
        }
    }

    return permTable;
}
let characterTableSpecht = memoised(characterTableSpechtHelper);


// Tensor product on class functions in the symmetric group on n letters.
function tensorProduct(n: number, left: Lin, right: Lin): Lin {
    let character = Lin.zero();
    for (let part of partitionsOf(n))
        character.set(part, left.get(part) * right.get(part));
    return character;
}

// Decompose a character into a linear combination of irreducible characters of S_n.
function decomposeCharacter(n: number, character: Lin): Lin {
    let irrCharacters = characterTableSpecht(n);
    return Lin.fromList(partitionsOf(n).map(partition => [partition, innerProduct(n, irrCharacters.get(partition), character)]))
}

// Return the trivial character for S_n.
function trivialCharacter(n: number): Lin {
    return Lin.fromList(partitionsOf(n).map(part => [part, 1]));
}

// Given a character 𝜒, return the character (g -> 𝜒(g^p)).
function characterPow(n: number, character: Lin, p: number): Lin {
    return Lin.fromList(partitionsOf(n).map(part => [part, character.get(cyclePow(part, p))]));
}

// Compute the (r + 1) exterior powers 𝛬^0(V), 𝛬^1(V), ..., 𝛬^r(V)
function exteriorPowers(n: number, r: number, character: Lin): Lin[] {
    let exteriors = [trivialCharacter(n)];

    // Now apply the recurrence e_i = (1/i)(e_{i-1} p_1 - e_{i-2} p_2 + ... + (-1)^(i-1) e_0 p_i)
    for (let i = 1; i <= r; i++) {
        let altSum = Lin.zero();
        for (let j = 1; j <= i; j++) {
            let term = tensorProduct(n, exteriors[i - j], characterPow(n, character, j));
            altSum = Lin.add(altSum, Lin.scale((j % 2 == 0) ? -1 : 1, term));
        }
        exteriors[i] = Lin.scale(1 / i, altSum);
    }

    return exteriors;
}

// Compute the (r + 1) symmetric powers S^0(V), S^1(V), ..., S^r(V)
function symmetricPowers(n: number, r: number, character: Lin): Lin[] {
    let symmetrics = [trivialCharacter(n)];

    // Now apply the recurrence e_i = (1/i)(e_{i-1} p_1 - e_{i-2} p_2 + ... + (-1)^(i-1) e_0 p_i)
    for (let i = 1; i <= r; i++) {
        let sum = Lin.zero();
        for (let j = 1; j <= i; j++) {
            let term = tensorProduct(n, symmetrics[i - j], characterPow(n, character, j));
            sum = Lin.add(sum, term);
        }
        symmetrics[i] = Lin.scale(1 / i, sum);
    }

    return symmetrics;
}

function tensorPowers(n: number, r: number, character: Lin): Lin[] {
    let tensors = [trivialCharacter(n)];
    for (let i = 1; i <= r; i++)
        tensors[i] = tensorProduct(n, tensors[i-1], character);
    return tensors;
}






// The rest of this file is just interfacing to HTML.







function showDomTable(n: number, table: CharacterTable, letter: string): HTMLElement {
    class Cl { constructor(readonly classname: string) {} }
    function C(classname: string) { return new Cl(classname); }
    let E = function(el: string, ...children: (string | Cl | HTMLElement)[]): HTMLElement {
        let $el = document.createElement(el);
        for (let $child of children) {
            if (typeof $child == 'string') {
                let $text = document.createTextNode($child);
                $el.appendChild($text);
            } else if ($child instanceof Cl)    {
                $el.classList.add($child.classname)
            } else{
                $el.appendChild($child);
            }
        }
        return $el;
    }
    let M = function(text: string): HTMLSpanElement {
        let $span = E('span');
        $span.innerHTML = katex.renderToString(text);
        return $span;
    }

    let selectedPartitions: string[] = [];

    let parts = partitionsOf(n);
    let partsReverse = parts.slice().reverse();

    let $table = E('table',
        E('thead',
            E('tr',
                C('underlined'),
                E('td', C('rightderlined')),
                ...partsReverse.map(part => E('th', M(`{${showPartitionCompact(part)}}`))))),
        E('tbody',
            // Row showing the sizes of conjugacy classes.
            E('tr',
                C('underlined'),
                E('td', C('rightderlined'), '#'),
                ...partsReverse.map(part => E('td', M("" + conjugacySize(part))))),

            // A row for each character. Each row is clickable, and will XOR itself onto the
            // current product of characters.
            ...parts.map(part => {
                let cells = partsReverse.map(colpart => E('td', M("" + table.get(part).get(colpart))));
                let row = E('tr',
                    E('td',
                        C('rightderlined'),
                        M(letter + showPartitionCompact(part))),
                    ...cells);
                row.classList.add('character');
                row.addEventListener('click', function() {
                    let partString = showPartition(part);
                    let idx = selectedPartitions.indexOf(partString);
                    if (idx >= 0) {
                        selectedPartitions.splice(idx, 1);
                        this.classList.remove('selected');
                    } else {
                        selectedPartitions.push(partString);
                        this.classList.add('selected');
                    }

                    updateProduct();
                });
                return row
            })));

    let $product = E('div');

    function productContents(): HTMLElement[] {
        if (selectedPartitions.length == 0)
            return [E('p', "Select some characters to tensor.")];

        let selectedCharacters = selectedPartitions.map(partitionString => table.get(readPartition(partitionString)));
        let chi = selectedCharacters.reduce((a, b) => tensorProduct(n, a, b));

        let interestingCharacters: [string, Lin][] = [
            ["\\chi", chi],
            ...exteriorPowers(n, 4, chi).slice(2).map((lin, r) => <[string, Lin]>[`\\wedge^{${r+2}} \\chi`, lin]),
            ...symmetricPowers(n, 4, chi).slice(2).map((lin, r) => <[string, Lin]>[`S^{${r+2}} \\chi`, lin]),
            ...tensorPowers(n, 4, chi).slice(2).map((lin, r) => <[string, Lin]>[`\\otimes^{${r+2}} \\chi`, lin])
        ];

        return [
            E('p',
                "Selected ",
                M(`\\chi = ${selectedPartitions.map(s => letter + s).join(` \\otimes `)}`)),
            E('h3',
                "Decomposition of ",
                M(`\\chi`),
                " and its exterior, symmetric and tensor powers"),
            createDecompositionTable(interestingCharacters)
        ]
    }

    function updateProduct() {
        while ($product.firstChild)
            $product.removeChild($product.firstChild);

        $product.append(...productContents());
    }

    function createDecompositionTable(characters: [string, Lin][]): HTMLElement {
        let irreducibles = characterTableSpecht(n);
        return E('table',
            // Header: the names of each character.
            E('tr',
                C('underlined'),
                E('td', C('rightderlined')),
                ...characters.map(([name, _]) => E('th', M(name)))),

            // Second row: dimensions.
            E('tr',
                C('underlined'),
                E('td', C('rightderlined'), M('\\dim')),
                ...characters.map(([_, lin]) => E('td', "" + lin.get(irreducibles.partitions[irreducibles.partitions.length - 1])))),

            // Rows: multiplicities.
            ...irreducibles.keys().map(partition => E('tr',
                E('td',
                    C('rightderlined'),
                    M(`s${showPartitionCompact(partition)}`)),
                ...characters.map(
                    ([_, lin]) => E('td', "" + innerProduct(n, irreducibles.get(partition), lin)))
            ))
        )
    }

    updateProduct();

    return E('div',
        E('h3',
            `Character table for the ${(letter == 's') ? 'Specht' : 'permutation'} modules of `,
            M(`S_{${n}}`)),
        $table,
        $product);
}

function doComputation() {
    let checked = <HTMLInputElement>document.querySelector('input[name="module-kind"]:checked');
    let n = (<HTMLInputElement>document.getElementById("order")!).valueAsNumber;
    let $tablePlace = <HTMLDivElement>document.getElementById("tablePlace");

    while ($tablePlace.firstChild)
        $tablePlace.removeChild($tablePlace.firstChild);

    let characterTable = (checked.value == "specht") ? characterTableSpecht(n) : characterTablePermutation(n);
    let letter = (checked.value == "specht") ? "s" : "h";

    $tablePlace.appendChild(showDomTable(n, characterTable, letter));
}

let $computationForm = <HTMLFormElement>document.getElementById('computationForm');
$computationForm.addEventListener('submit', ev => ev.preventDefault());
$computationForm.addEventListener('change', ev => {ev.preventDefault(); doComputation();});
doComputation();
