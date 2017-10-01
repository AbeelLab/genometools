package abeel.genometools.gfa

import scala.collection.mutable.MutableList

case class Segment(val idx: Int, val sequenceLen: Int, val genomeIdx: List[Int], val incoming: MutableList[Int], val outgoing: MutableList[Int]) {
  override def toString() = {
    idx + "," + sequenceLen + "," + genomeIdx + "," + incoming + "," + outgoing
  }
}