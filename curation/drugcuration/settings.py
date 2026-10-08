"""Scrapy settings: polite, cached, resumable. Per-spider delays live in each spider's custom_settings."""
import os

BOT_NAME = "drugcuration"
SPIDER_MODULES = ["drugcuration.spiders"]
NEWSPIDER_MODULE = "drugcuration.spiders"

# Identify honestly. Put your own contact address in the env var CURATION_CONTACT.
USER_AGENT = f"DrugCurationBot/0.1 (academic non-commercial research; contact: {os.environ.get('CURATION_CONTACT', 'unset')})"
ROBOTSTXT_OBEY = True

CONCURRENT_REQUESTS = 4
CONCURRENT_REQUESTS_PER_DOMAIN = 1
DOWNLOAD_DELAY = 1.0
AUTOTHROTTLE_ENABLED = True
AUTOTHROTTLE_START_DELAY = 1.0
AUTOTHROTTLE_MAX_DELAY = 30.0
AUTOTHROTTLE_TARGET_CONCURRENCY = 1.0

RETRY_ENABLED = True
RETRY_TIMES = 3
RETRY_HTTP_CODES = [429, 500, 502, 503, 504, 522, 524, 408]

# Every page is cached on disk, so re-runs never hit the sites again.
HTTPCACHE_ENABLED = True
HTTPCACHE_DIR = "httpcache"
HTTPCACHE_EXPIRATION_SECS = 0
HTTPCACHE_IGNORE_HTTP_CODES = [403, 429, 500, 502, 503, 504]

COOKIES_ENABLED = True
TELNETCONSOLE_ENABLED = False
LOG_LEVEL = "INFO"
FEED_EXPORT_ENCODING = "utf-8"
REQUEST_FINGERPRINTER_IMPLEMENTATION = "2.7"
TWISTED_REACTOR = "twisted.internet.asyncioreactor.AsyncioSelectorReactor"

ITEM_PIPELINES = {"drugcuration.pipelines.DomainLogPipeline": 100}

# ---------------------------------------------------------------------------
# scrapy-playwright (https://github.com/scrapy-plugins/scrapy-playwright)
# go.drugbank.com answers plain HTTP clients with a JavaScript "Securing Connection" challenge,
# so the DrugBank spider renders pages in headless Chromium. Only requests whose meta contains
# {"playwright": True} go through the browser; ChEMBL / openFDA / WHOCC stay on plain Scrapy HTTP.
# Set USE_PLAYWRIGHT=0 to switch the browser off completely.
# ---------------------------------------------------------------------------
USE_PLAYWRIGHT = os.environ.get("USE_PLAYWRIGHT", "1") != "0"
if USE_PLAYWRIGHT:
    DOWNLOAD_HANDLERS = {
        "http": "scrapy_playwright.handler.ScrapyPlaywrightDownloadHandler",
        "https": "scrapy_playwright.handler.ScrapyPlaywrightDownloadHandler",
    }
    PLAYWRIGHT_BROWSER_TYPE = "chromium"
    # PLAYWRIGHT_HEADLESS=0 opens a visible browser window (useful if a challenge needs a human click).
    PLAYWRIGHT_LAUNCH_OPTIONS = {"headless": os.environ.get("PLAYWRIGHT_HEADLESS", "1") != "0"}
    PLAYWRIGHT_DEFAULT_NAVIGATION_TIMEOUT = 60_000  # ms
    PLAYWRIGHT_MAX_PAGES_PER_CONTEXT = 2
    # Persistent browser profile: cookies (including a passed site check) survive between runs, exactly
    # like a normal browser. Delete the .pw-profile folder to start fresh.
    PLAYWRIGHT_CONTEXTS = {"default": {"user_data_dir": os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ".pw-profile")}}
    PLAYWRIGHT_MAX_CONTEXTS = 1

    def _abort_heavy(request):
        """Skip images, fonts, media and trackers: faster pages, less load on the site."""
        return request.resource_type in {"image", "media", "font"} or any(
            t in request.url for t in ("google-analytics", "googletagmanager", "doubleclick", "hotjar"))

    PLAYWRIGHT_ABORT_REQUEST = _abort_heavy
