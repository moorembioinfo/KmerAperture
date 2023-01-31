(*
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*)

(* (c) 2023 Paolo Ribeca, <paolo.ribeca@gmail.com>
    Compile with
     ocamlopt.opt -O3 -o KmerApertureParser KmerApertureParser.ml -ccopt -static
    Parameters are file name and k-mer size *)

let rc s =
  let b = Bytes.of_string s in
  let compl = function
    | 'A' | 'a' -> 'T'
    | 'C' | 'c' -> 'G'
    | 'G' | 'g' -> 'C'
    | 'T' | 't' -> 'A'
    | _ -> 'N' in
  let len = Bytes.length b in
  let red_len = len - 1 in
  let half_len = red_len / 2 in
  for i = 0 to half_len do
    let idx = red_len - i in
    let c = Bytes.get b i in
    Bytes.set b i (compl (Bytes.get b idx));
    Bytes.set b idx (compl c)
  done;
  Bytes.to_string b

let () =
  if Array.length Sys.argv <> 3 then begin
    Printf.eprintf "Syntax: KmerApertureParser <FASTA_input> <k>\n";
    Printf.eprintf " The input file can contain more than one sequence\n%!";
    exit 1
  end;
  let file = open_in Sys.argv.(1) and k = int_of_string Sys.argv.(2) in
  let buf = Buffer.create 1024 in
  let process_contig () =
    let seq = Buffer.contents buf in
    let top = String.length seq - k in
    for i = 0 to top do
      let kmer = String.sub seq i k in
      let kmer_rc = rc kmer in
      min kmer kmer_rc |> Printf.printf "%s\n"
    done;
    Buffer.clear buf in
  begin try
    while true do
      let line = input_line file in
      if line <> "" then begin
        if line.[0] = '>' then
          process_contig ()
        else
          Buffer.add_string buf line
      end
    done
  with End_of_file ->
    process_contig ()
  end;
  close_in file

