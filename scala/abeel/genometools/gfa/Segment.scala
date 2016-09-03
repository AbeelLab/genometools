package abeel.genometools.gfa

import scala.collection.mutable.MutableList

case class Segment(val idx: Int, val sequence: String, val genomes: List[String], val incoming: MutableList[Int], val outgoing: MutableList[Int]) {
  override def toString() = {
    idx + "," + (if (sequence.size > 3) sequence.take(3) + "..." else sequence) + "," + genomes + "," + incoming + "," + outgoing
  }
}