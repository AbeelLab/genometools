package abeel.genometools.gfa

import scala.collection.mutable.HashMap

object Dictionary {

  private val map = new HashMap[String, Int];

  def get(key: String) = {
    if (!map.contains(key))
      map += key -> map.size
    map(key)
  }

  def keyMap = map

}