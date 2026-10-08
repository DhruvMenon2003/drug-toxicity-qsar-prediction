"""Adds the domain of every URL an item came from to `domains_visited` (feeds the Filtered_Domains_Log column)."""
from urllib.parse import urlparse


class DomainLogPipeline:
    def process_item(self, item, *args):  # Scrapy <2.14 passes spider; newer does not
        urls = item.get("source_urls") or []
        item["domains_visited"] = sorted({urlparse(u).netloc for u in urls if u})
        return item
